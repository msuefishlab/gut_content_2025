#!/usr/bin/env bash
# =============================================================================
# parasite_rrna_annotation.sh
# rDNA repeat unit annotation for assembled rRNA locus contigs
# =============================================================================
# Purpose: Annotate an assembled rRNA locus contig by:
#   1. Running Barrnap to locate rRNA features and identify repeat unit boundaries
#   2. Splitting the contig into individual rDNA repeat units
#   3. Reverse-complementing minus-strand units so ITSx sees canonical orientation
#   4. Running ITSx on each unit to extract SSU/ITS1/5.8S/ITS2/LSU regions
#   5. BLASTn on extracted SSU and ITS regions for species identification
#
# Background: rDNA exists as a tandem repeat array (18S-ITS1-5.8S-ITS2-28S).
#   Assembled contigs often span multiple units. Feeding a multi-unit contig
#   to ITSx causes it to report inflated ITS lengths spanning unit boundaries.
#   This script uses Barrnap coordinates to split before running ITSx.
#
# Usage:
#   bash parasite_rrna_annotation.sh -i ASSEMBLY.fa -c CONTIG_ID [options]
#
# Required:
#   -i  PATH    Assembly FASTA (may contain multiple contigs)
#   -c  STR     Contig ID to annotate (e.g. contig_11)
#
# Optional:
#   -o  PATH    Output directory [default: rrna_annotation]
#   -e  STR     Conda environment containing barrnap and ITSx [default: mitos2]
#   -k  STR     Barrnap kingdom: euk|bac|mito [default: euk]
#   -n  PATH    Local BLAST nr database (omit = remote NCBI)
#   -E  FLOAT   BLASTn e-value [default: 1e-10]
#   -t  INT     Max BLASTn targets [default: 10]
#   -T  INT     Threads [default: 4]
#   -h          Show this help
#
# Dependencies:
#   samtools, python3, conda (with barrnap + ITSx in -e environment), blastn
#
# Output:
#   $OUTDIR/barrnap.gff              — raw barrnap annotations
#   $OUTDIR/units/unit_N.fa          — individual rDNA repeat unit FASTAs
#   $OUTDIR/itsx/unit_N.*            — ITSx output per unit
#   $OUTDIR/itsx/unit_N.positions.txt — ITSx region coordinates
#   $OUTDIR/blast/SSU_blast.tsv      — BLASTn hits for 18S sequences
#   $OUTDIR/blast/ITS_blast.tsv      — BLASTn hits for ITS1+5.8S+ITS2
#   $OUTDIR/summary.txt              — region sizes and top BLAST hit per unit
# =============================================================================

set -euo pipefail

# =============================================================================
# Defaults and argument parsing
# =============================================================================

ASSEMBLY=""
CONTIG_ID=""
OUTDIR="rrna_annotation"
CONDA_ENV="mitos2"
KINGDOM="euk"
NR_DB=""
EVALUE="1e-10"
MAX_SEQS=10
THREADS=4

usage() {
    sed -n '16,42p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 1
}

while getopts "i:c:o:e:k:n:E:t:T:h" opt; do
    case $opt in
        i) ASSEMBLY="$OPTARG" ;;
        c) CONTIG_ID="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        k) KINGDOM="$OPTARG" ;;
        n) NR_DB="$OPTARG" ;;
        E) EVALUE="$OPTARG" ;;
        t) MAX_SEQS="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$ASSEMBLY" || -z "$CONTIG_ID" ]]; then
    echo "ERROR: -i (assembly) and -c (contig ID) are required."
    usage
fi

mkdir -p "${OUTDIR}/units" "${OUTDIR}/itsx" "${OUTDIR}/blast"
LOG="${OUTDIR}/pipeline.log"
echo "[$(date)] rRNA annotation started" | tee "$LOG"
echo "[$(date)] Assembly: $ASSEMBLY  Contig: $CONTIG_ID" | tee -a "$LOG"

# =============================================================================
# STEP 1: Extract target contig and run Barrnap
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 1: Extract contig and run Barrnap ===" | tee -a "$LOG"

CONTIG_FA="${OUTDIR}/${CONTIG_ID}.fa"
samtools faidx "$ASSEMBLY" "$CONTIG_ID" > "$CONTIG_FA"
contig_len=$(grep -v "^>" "$CONTIG_FA" | tr -d '\n' | wc -c | tr -d ' ')
echo "[$(date)] Contig length: ${contig_len} bp" | tee -a "$LOG"

GFF="${OUTDIR}/barrnap.gff"
conda run -n "$CONDA_ENV" barrnap --kingdom "$KINGDOM" "$CONTIG_FA" \
    > "$GFF" 2>>"$LOG"

echo "[$(date)] Barrnap features:" | tee -a "$LOG"
grep -v "^#" "$GFF" | awk '{printf "  %s-%s  %s  strand=%s  score=%s\n", $4,$5,$9,$7,$6}' \
    | tee -a "$LOG"

# =============================================================================
# STEP 2: Parse Barrnap GFF and determine repeat unit boundaries
# =============================================================================
# Strategy: each rDNA repeat unit is anchored by an SSU (18S). When multiple
# SSU features appear on the contig, each one marks the boundary of a unit.
# For minus-strand units the SSU appears after the LSU in coordinate space,
# so we sort all SSU positions and use them as split points.
#
# Units that lack a detected SSU (e.g. junction fragments from mid-contig
# splits) are kept as-is but flagged as "partial" in the summary.

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 2: Determine repeat unit boundaries ===" | tee -a "$LOG"

python3 - <<PYEOF
import sys, os

gff_file   = "${GFF}"
contig_id  = "${CONTIG_ID}"
contig_len = int("${contig_len}")
outdir     = "${OUTDIR}"
contig_fa  = "${CONTIG_FA}"

features = []
with open(gff_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        seq_id, _, feat_type, start, end, score, strand, _, attrs = parts
        name = [a.split("=")[1] for a in attrs.split(";") if a.startswith("Name=")]
        name = name[0] if name else feat_type
        features.append({
            "name": name, "start": int(start), "end": int(end),
            "strand": strand, "score": score
        })

# Find SSU features — they anchor each repeat unit
ssu_feats = [f for f in features if "18S_rRNA" in f["name"]]

if len(ssu_feats) == 0:
    # No SSU found: treat entire contig as one unit
    boundaries = [(1, contig_len)]
    print("  No SSU features found — treating contig as single unit")
elif len(ssu_feats) == 1:
    boundaries = [(1, contig_len)]
    print(f"  Single SSU at {ssu_feats[0]['start']}-{ssu_feats[0]['end']} "
          f"({ssu_feats[0]['strand']}) — single unit")
else:
    # Multiple SSUs: determine boundaries between them.
    # For a minus-strand unit: order is LSU...5.8S...SSU (high→low coords),
    # so the SSU end coordinate is the highest position in the unit.
    # For a plus-strand unit: SSU start is the lowest position.
    # Simple approach: collect the midpoints between consecutive SSU features
    # (sorted by start coordinate) as split points.
    ssu_sorted = sorted(ssu_feats, key=lambda f: f["start"])
    splits = [1]
    for i in range(1, len(ssu_sorted)):
        prev_end   = ssu_sorted[i-1]["end"]
        curr_start = ssu_sorted[i]["start"]
        # Split in the gap between the end of one SSU and the start of the next
        midpoint = (prev_end + curr_start) // 2
        splits.append(midpoint)
    splits.append(contig_len)
    boundaries = [(splits[i] + (1 if i > 0 else 0), splits[i+1])
                  for i in range(len(splits)-1)]
    print(f"  Found {len(ssu_feats)} SSU features — splitting into {len(boundaries)} units")
    for i, (s, e) in enumerate(boundaries, 1):
        print(f"    unit{i}: positions {s}-{e} ({e-s+1} bp)")

# Write boundaries file for shell to read
bounds_file = os.path.join(outdir, "unit_boundaries.txt")
with open(bounds_file, "w") as fh:
    for i, (s, e) in enumerate(boundaries, 1):
        fh.write(f"unit{i}\t{s}\t{e}\n")

print(f"  Boundaries written to {bounds_file}")
PYEOF

echo "" | tee -a "$LOG"

# =============================================================================
# STEP 3: Extract units and reverse-complement minus-strand units
# =============================================================================
# ITSx expects the canonical 18S→28S orientation (plus strand).
# If Barrnap's strongest SSU hit for a unit is on the minus strand,
# we reverse-complement the unit so ITSx can parse it correctly.

echo "[$(date)] === STEP 3: Extract and orient repeat units ===" | tee -a "$LOG"

while IFS=$'\t' read -r unit_name start end; do
    out_fa="${OUTDIR}/units/${unit_name}.fa"
    samtools faidx "$CONTIG_FA" "${CONTIG_ID}:${start}-${end}" \
        | sed "s/>.*/>${unit_name}/" > "$out_fa"
    unit_len=$(grep -v "^>" "$out_fa" | tr -d '\n' | wc -c | tr -d ' ')

    # Determine dominant strand from barrnap SSU features in this region
    dominant_strand=$(python3 - <<PYEOF2
import sys
features = []
with open("${GFF}") as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9: continue
        name = [a.split("=")[1] for a in parts[8].split(";") if a.startswith("Name=")]
        name = name[0] if name else ""
        if "18S_rRNA" in name:
            s, e = int(parts[3]), int(parts[4])
            if s >= $start and e <= $end:
                features.append(parts[6])
# majority strand of SSU features in this window
if not features:
    print("+")
elif features.count("-") > features.count("+"):
    print("-")
else:
    print("+")
PYEOF2
)

    if [[ "$dominant_strand" == "-" ]]; then
        # Reverse complement using python (avoids seqtk dependency)
        python3 - <<PYEOF3
seq = open("${out_fa}").read().split("\n")
header = seq[0]
dna = "".join(seq[1:])
comp = str.maketrans("ACGTacgt", "TGCAtgca")
rc = dna.translate(comp)[::-1]
# rewrap at 60 chars
wrapped = "\n".join(rc[i:i+60] for i in range(0, len(rc), 60))
with open("${out_fa}", "w") as fh:
    fh.write(header + "_rc\n" + wrapped + "\n")
PYEOF3
        echo "[$(date)] ${unit_name}: ${unit_len}bp (minus strand → reverse-complemented)" | tee -a "$LOG"
    else
        echo "[$(date)] ${unit_name}: ${unit_len}bp (plus strand)" | tee -a "$LOG"
    fi

done < "${OUTDIR}/unit_boundaries.txt"

# =============================================================================
# STEP 4: Run ITSx on each unit
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 4: ITSx annotation per unit ===" | tee -a "$LOG"

while IFS=$'\t' read -r unit_name start end; do
    unit_fa="${OUTDIR}/units/${unit_name}.fa"
    itsx_prefix="${OUTDIR}/itsx/${unit_name}"

    echo "[$(date)] ITSx: ${unit_name}..." | tee -a "$LOG"
    conda run -n "$CONDA_ENV" ITSx \
        -i "$unit_fa" \
        -o "$itsx_prefix" \
        --taxa Metazoa \
        --save_regions all \
        --complement T \
        --positions T \
        --cpu "$THREADS" \
        2>>"$LOG" || true

    pos_file="${itsx_prefix}.positions.txt"
    if [[ -f "$pos_file" ]]; then
        cat "$pos_file" | tee -a "$LOG"
    else
        echo "[$(date)]   ${unit_name}: no ITSx positions output" | tee -a "$LOG"
    fi

done < "${OUTDIR}/unit_boundaries.txt"

# =============================================================================
# STEP 5: BLASTn SSU and ITS regions for species identification
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 5: BLASTn species identification ===" | tee -a "$LOG"

BLAST_FMT="6 qseqid sseqid pident length qcovs evalue bitscore stitle sscinames"

# Collect SSU and ITS sequences across all units into combined query files
SSU_FA="${OUTDIR}/blast/all_SSU.fa"
ITS_FA="${OUTDIR}/blast/all_ITS.fa"
: > "$SSU_FA"
: > "$ITS_FA"

while IFS=$'\t' read -r unit_name start end; do
    itsx_prefix="${OUTDIR}/itsx/${unit_name}"
    ssu_f="${itsx_prefix}.SSU.fasta"
    its1_f="${itsx_prefix}.ITS1.fasta"
    s58_f="${itsx_prefix}.5_8S.fasta"
    its2_f="${itsx_prefix}.ITS2.fasta"

    [[ -s "$ssu_f" ]] && cat "$ssu_f" >> "$SSU_FA"

    # Concatenate ITS1 + 5.8S + ITS2 into a single ITS region per unit
    if [[ -s "$its1_f" && -s "$s58_f" && -s "$its2_f" ]]; then
        # Merge the three regions into one sequence using python
        python3 - <<PYEOF4
import re
def read_fa(path):
    lines = open(path).read().strip().split("\n")
    return "".join(lines[1:])
its1 = read_fa("$its1_f")
s58  = read_fa("$s58_f")
its2 = read_fa("$its2_f")
merged = its1 + s58 + its2
wrapped = "\n".join(merged[i:i+60] for i in range(0, len(merged), 60))
with open("$ITS_FA", "a") as fh:
    fh.write(f">${unit_name}_ITS_region\n{wrapped}\n")
PYEOF4
    fi

done < "${OUTDIR}/unit_boundaries.txt"

# Run BLASTn on whichever query files have sequences
for query_fa in "$SSU_FA" "$ITS_FA"; do
    [[ -s "$query_fa" ]] || continue
    label=$(basename "$query_fa" .fa)
    out_tsv="${OUTDIR}/blast/${label}.blast.tsv"
    out_txt="${OUTDIR}/blast/${label}.blast.txt"

    if [[ -s "$out_tsv" ]]; then
        echo "[$(date)] BLASTn: ${label} — already done, skipping." | tee -a "$LOG"
    else
        echo "[$(date)] BLASTn: ${label}..." | tee -a "$LOG"
        if [[ -n "$NR_DB" ]]; then
            blastn -query "$query_fa" -db "$NR_DB" \
                -out "$out_tsv" -outfmt "$BLAST_FMT" \
                -evalue "$EVALUE" -max_target_seqs "$MAX_SEQS" \
                -num_threads "$THREADS" 2>>"$LOG"
            blastn -query "$query_fa" -db "$NR_DB" \
                -out "$out_txt" -outfmt 0 \
                -evalue "$EVALUE" -max_target_seqs 3 \
                -num_threads "$THREADS" 2>>"$LOG"
        else
            echo "[$(date)]   NOTE: using remote BLAST — may take several minutes." | tee -a "$LOG"
            blastn -query "$query_fa" -db nr \
                -out "$out_tsv" -outfmt "$BLAST_FMT" \
                -evalue "$EVALUE" -max_target_seqs "$MAX_SEQS" \
                -remote 2>>"$LOG"
            blastn -query "$query_fa" -db nr \
                -out "$out_txt" -outfmt 0 \
                -evalue "$EVALUE" -max_target_seqs 3 \
                -remote 2>>"$LOG"
        fi
    fi

    if [[ -s "$out_tsv" ]]; then
        echo "[$(date)]   Top hits:" | tee -a "$LOG"
        head -3 "$out_tsv" | awk -F'\t' \
            '{printf "  %s -> %s | pident=%.1f%% qcov=%s%% evalue=%s\n",
              $1, $9, $3, $5, $6}' | tee -a "$LOG"
    else
        echo "[$(date)]   No hits above evalue ${EVALUE}." | tee -a "$LOG"
    fi
done

# =============================================================================
# STEP 6: Summary
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 6: Summary ===" | tee -a "$LOG"

SUMMARY="${OUTDIR}/summary.txt"
{
printf "%-10s %-8s %-10s %-10s %-10s %-10s %-10s  %s\n" \
    "Unit" "Length" "SSU" "ITS1" "5.8S" "ITS2" "LSU" "Top BLAST hit"
printf "%-10s %-8s %-10s %-10s %-10s %-10s %-10s  %s\n" \
    "----" "------" "---" "----" "----" "----" "---" "-------------"

while IFS=$'\t' read -r unit_name start end; do
    unit_len=$(( end - start + 1 ))
    pos_file="${OUTDIR}/itsx/${unit_name}.positions.txt"

    # Parse ITSx positions — use grep/sed (BSD awk lacks 3-arg match())
    ssu="" its1="" s58="" its2="" lsu=""
    if [[ -f "$pos_file" ]]; then
        ssu=$(grep -oE  'SSU: [0-9]+-[0-9]+|SSU: Not found'  "$pos_file" | sed 's/SSU: //'   || true)
        its1=$(grep -oE 'ITS1: [0-9]+-[0-9]+|ITS1: Not found' "$pos_file" | sed 's/ITS1: //'  || true)
        s58=$(grep -oE  '5\.8S: [0-9]+-[0-9]+|5\.8S: Not found' "$pos_file" | sed 's/5\.8S: //' || true)
        its2=$(grep -oE 'ITS2: [0-9]+-[0-9]+|ITS2: Not found' "$pos_file" | sed 's/ITS2: //'  || true)
        lsu=$(grep -oE  'LSU: [0-9]+-[0-9]+|LSU: Not found'   "$pos_file" | sed 's/LSU: //'   || true)
    fi

    # Top BLAST hit for this unit (SSU preferred, then ITS)
    top_hit="—"
    for tsv in "${OUTDIR}/blast/all_SSU.blast.tsv" "${OUTDIR}/blast/all_ITS.blast.tsv"; do
        [[ -f "$tsv" ]] || continue
        h=$(grep "^${unit_name}" "$tsv" | head -1 | awk -F'\t' \
            '{printf "%s (%.1f%%)", $9, $3}' 2>/dev/null || true)
        [[ -n "$h" ]] && top_hit="$h" && break
    done

    printf "%-10s %-8s %-10s %-10s %-10s %-10s %-10s  %s\n" \
        "$unit_name" "${unit_len}bp" \
        "${ssu:-N/A}" "${its1:-N/A}" "${s58:-N/A}" \
        "${its2:-N/A}" "${lsu:-N/A}" "$top_hit"

done < "${OUTDIR}/unit_boundaries.txt"
} | tee "$SUMMARY" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "[$(date)] Done." | tee -a "$LOG"
echo "Summary:    $SUMMARY"
echo "ITSx files: ${OUTDIR}/itsx/"
echo "BLAST hits: ${OUTDIR}/blast/"
