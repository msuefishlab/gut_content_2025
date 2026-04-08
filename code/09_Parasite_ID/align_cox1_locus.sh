#!/usr/bin/env bash
# align_cox1_locus.sh
# Align recruited ONT reads and Eustrongylides reference sequences to the
# assembled COX1 locus (mtDNA_02 / contig_3) to diagnose the 13 bp gap
# region (contig_3:2962-2974) between MITOS2-annotated exons.
#
# Outputs BAM files for IGV inspection and a coverage report at the gap.
#
# Usage: bash align_cox1_locus.sh [options]
#   -r  PATH    Recruited ONT reads FASTA
#   -a  PATH    Assembly FASTA (mtDNA_02)
#   -g  PATH    Reference sequences GenBank file
#   -c  STR     Contig ID to focus on [default: contig_3]
#   -o  PATH    Output directory
#   -T  INT     Threads [default: 4]
#   -h          Help

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve repository root (script lives in code/09_Parasite_ID/)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
READS=""
ASSEMBLY=""
REFGB="${ROOT}/input_data/09_Parasite_ID/eustrongylides_cox1_refs.gb"
CONTIG_ID="contig_3"
OUTDIR=""
THREADS=4

GAP_START=2962
GAP_END=2974
CONTEXT=50   # bp of flanking context for depth report

usage() {
    cat <<EOF
Usage: bash $(basename "$0") [options]

Aligns recruited ONT reads and Eustrongylides COX1 reference sequences to the
assembled COX1 locus to diagnose the gap at ${CONTIG_ID}:${GAP_START}-${GAP_END}.

Options:
  -r  PATH    Recruited ONT reads FASTA
                [default: auto-detect reads_by_target/reads_mtDNA_02.fa]
  -a  PATH    Assembly FASTA (mtDNA_02)
                [default: auto-detect assemblies/mtDNA_02/assembly.fasta]
  -g  PATH    Reference sequences GenBank file
                [default: ${REFGB}]
  -c  STR     Contig ID to focus on [default: ${CONTIG_ID}]
  -o  PATH    Output directory
                [default: output_data/09_Parasite_ID/<run>/cox1_alignment]
  -T  INT     Threads [default: ${THREADS}]
  -h          Show this help and exit
EOF
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts "r:a:g:c:o:T:h" opt; do
    case "$opt" in
        r) READS="$OPTARG" ;;
        a) ASSEMBLY="$OPTARG" ;;
        g) REFGB="$OPTARG" ;;
        c) CONTIG_ID="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Auto-detect paths if not provided
# ---------------------------------------------------------------------------
detect_run_dir() {
    local base="${ROOT}/output_data/09_Parasite_ID"
    # Find most recent run directory (april*_run or similar)
    local run
    run=$(find "$base" -maxdepth 1 -type d -name '*_run' \
        | sort | tail -1)
    if [[ -z "$run" ]]; then
        run=$(find "$base" -maxdepth 1 -type d | grep -v "^${base}$" | sort | tail -1)
    fi
    echo "$run"
}

if [[ -z "$READS" || -z "$ASSEMBLY" || -z "$OUTDIR" ]]; then
    RUN_DIR=$(detect_run_dir)
    if [[ -z "$RUN_DIR" ]]; then
        echo "ERROR: Could not auto-detect run directory under output_data/09_Parasite_ID/" >&2
        exit 1
    fi
fi

[[ -z "$READS" ]]    && READS="${RUN_DIR}/reads_by_target/reads_mtDNA_02.fa"
[[ -z "$ASSEMBLY" ]] && ASSEMBLY="${RUN_DIR}/assemblies/mtDNA_02/assembly.fasta"
[[ -z "$OUTDIR" ]]   && OUTDIR="${RUN_DIR}/cox1_alignment"

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
echo "================================================================="
echo "  COX1 Locus Alignment Diagnostic"
echo "  $(date)"
echo "================================================================="
echo "  Reads:    $READS"
echo "  Assembly: $ASSEMBLY"
echo "  Refs:     $REFGB"
echo "  Contig:   $CONTIG_ID"
echo "  Output:   $OUTDIR"
echo "  Threads:  $THREADS"
echo "================================================================="
echo ""

for f in "$READS" "$ASSEMBLY" "$REFGB"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
        exit 1
    fi
done

for tool in minimap2 samtools python3 blastn; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: Required tool not found in PATH: $tool" >&2
        exit 1
    fi
done

python3 -c "from Bio import SeqIO" 2>/dev/null \
    || { echo "ERROR: Biopython not available. Install with: pip install biopython" >&2; exit 1; }

mkdir -p "${OUTDIR}/refs" "${OUTDIR}/alignments"

REPORT="${OUTDIR}/cox1_alignment_report.txt"
LOG="${OUTDIR}/align_cox1.log"

{
    echo "COX1 Locus Alignment Report"
    echo "Generated: $(date)"
    echo "Assembly:  $ASSEMBLY"
    echo "Reads:     $READS ($(grep -c '^>' "$READS") sequences)"
    echo "Refs:      $REFGB"
    echo ""
} > "$REPORT"

echo "[$(date)] Starting COX1 alignment pipeline" | tee -a "$LOG"

# ---------------------------------------------------------------------------
# STEP 1: Extract reference sequences from GenBank
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 1: Extract Eustrongylides reference sequences from GenBank"
echo "================================================================="

REFS_FA="${OUTDIR}/refs/eustrongylides_cox1_refs.fa"

python3 - <<PYEOF 2>&1 | tee -a "$LOG"
from Bio import SeqIO
import sys

gb_file  = "${REFGB}"
out_fa   = "${REFS_FA}"
report   = "${REPORT}"

records = list(SeqIO.parse(gb_file, "genbank"))
print(f"Parsed {len(records)} records from GenBank file")

header = f"{'Accession':<12}  {'Organism':<40}  {'Len':>6}  {'Country'}"
print(header)
print("-" * len(header))

with open(out_fa, "w") as fh, open(report, "a") as rh:
    rh.write("=== Reference Sequences ===\n")
    rh.write(header + "\n")
    rh.write("-" * len(header) + "\n")
    for rec in records:
        acc      = rec.id
        org      = rec.annotations.get("organism", "unknown")
        length   = len(rec.seq)
        country  = ""
        for feat in rec.features:
            if feat.type == "source":
                country = feat.qualifiers.get("country", [""])[0]
                break
        row = f"{acc:<12}  {org:<40}  {length:>6}  {country}"
        print(row)
        rh.write(row + "\n")
        fh.write(f">{acc} {org} ({length} bp)\n{str(rec.seq)}\n")
    rh.write(f"\nTotal: {len(records)} sequences written to {out_fa}\n\n")

print(f"\nWrote {len(records)} sequences to {out_fa}")
PYEOF

echo "[$(date)] STEP 1 done" | tee -a "$LOG"

# ---------------------------------------------------------------------------
# STEP 2: Align ONT reads to assembly
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 2: Align recruited ONT reads to assembly (minimap2 map-ont)"
echo "================================================================="

ONT_BAM="${OUTDIR}/alignments/cox1_ont_reads.bam"

if [[ -s "$ONT_BAM" && -s "${ONT_BAM}.bai" ]]; then
    echo "[$(date)] ONT BAM already exists — skipping alignment" | tee -a "$LOG"
else
    echo "[$(date)] Running minimap2 (map-ont) ..." | tee -a "$LOG"
    minimap2 -ax map-ont -t "$THREADS" \
        --secondary=no \
        "$ASSEMBLY" "$READS" \
        2>>"$LOG" \
    | samtools sort -@ "$THREADS" -o "$ONT_BAM"
    samtools index "$ONT_BAM"
    echo "[$(date)] ONT BAM written: $ONT_BAM" | tee -a "$LOG"
fi

# Alignment stats
echo "[$(date)] ONT alignment stats:" | tee -a "$LOG"
samtools flagstat "$ONT_BAM" | tee -a "$LOG"
{
    echo "=== ONT Read Alignment Stats ==="
    samtools flagstat "$ONT_BAM"
    echo ""
} >> "$REPORT"

# ---------------------------------------------------------------------------
# STEP 3: Align reference sequences to assembly (blastn → SAM/BAM)
# ---------------------------------------------------------------------------
# minimap2 fails at ~80% identity. Use blastn which handles this divergence,
# then convert the pairwise alignment to SAM format for IGV visualization.
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 3: Align Eustrongylides references to assembly (blastn → BAM)"
echo "================================================================="

REFS_BAM="${OUTDIR}/alignments/cox1_refs.bam"
BLAST_TSV="${OUTDIR}/alignments/refs_vs_assembly.blastn.tsv"

if [[ -s "$REFS_BAM" && -s "${REFS_BAM}.bai" ]]; then
    echo "[$(date)] Refs BAM already exists — skipping alignment" | tee -a "$LOG"
else
    # Run blastn with aligned sequences in output
    echo "[$(date)] Running blastn (task=blastn, perc_identity>=70) ..." | tee -a "$LOG"
    blastn \
        -query "$REFS_FA" \
        -subject "$ASSEMBLY" \
        -outfmt "6 qseqid sseqid qstart qend qlen sstart send sstrand qseq sseq pident" \
        -task blastn \
        -perc_identity 70 \
        -max_hsps 1 \
        -max_target_seqs 1 \
        2>>"$LOG" \
    > "$BLAST_TSV"

    n_hits=$(wc -l < "$BLAST_TSV")
    echo "[$(date)] blastn: $n_hits alignments" | tee -a "$LOG"

    # Convert blastn tabular output to BAM via Python
    REFS_SAM="${OUTDIR}/alignments/cox1_refs.sam"
    python3 - <<PYEOF 2>>"$LOG"
import subprocess, sys

blast_tsv = "${BLAST_TSV}"
assembly  = "${ASSEMBLY}"
refs_fa   = "${REFS_FA}"
out_sam   = "${REFS_SAM}"

# ── helpers ──────────────────────────────────────────────────────────────────

def pairwise_to_cigar(qseq_aln, sseq_aln, qstart, qend, qlen):
    """
    Convert a BLAST pairwise alignment to a SAM CIGAR string.
    qseq_aln / sseq_aln are the aligned strings (may contain '-').
    Unaligned ends of the query become soft clips.
    """
    ops = []
    for q, s in zip(qseq_aln, sseq_aln):
        if q == '-':
            ops.append('D')          # deletion in read relative to ref
        elif s == '-':
            ops.append('I')          # insertion in read relative to ref
        else:
            ops.append('M')
    # Build run-length CIGAR
    cigar = ''
    if not ops:
        return '*'
    prev, count = ops[0], 1
    for op in ops[1:]:
        if op == prev:
            count += 1
        else:
            cigar += f"{count}{prev}"
            prev, count = op, 1
    cigar += f"{count}{prev}"
    # Add soft clips for unaligned query ends (qstart/qend are 1-based)
    left_clip  = qstart - 1
    right_clip = qlen - qend
    if left_clip:
        cigar = f"{left_clip}S" + cigar
    if right_clip:
        cigar = cigar + f"{right_clip}S"
    return cigar

def rc(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

# ── load reference sequences ──────────────────────────────────────────────────
from Bio import SeqIO
ref_seqs = {}
for rec in SeqIO.parse(refs_fa, 'fasta'):
    ref_seqs[rec.id] = str(rec.seq)

# ── load assembly contig lengths for SAM header ───────────────────────────────
contig_lens = {}
for rec in SeqIO.parse(assembly, 'fasta'):
    contig_lens[rec.id] = len(rec.seq)

# ── build SAM ─────────────────────────────────────────────────────────────────
written = 0
with open(out_sam, 'w') as fh:
    # Header
    fh.write('@HD\tVN:1.6\tSO:unsorted\n')
    for cid, clen in contig_lens.items():
        fh.write(f'@SQ\tSN:{cid}\tLN:{clen}\n')
    fh.write('@PG\tID:blastn\tPN:blastn\n')

    with open(blast_tsv) as bf:
        for line in bf:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t')
            qseqid, sseqid = parts[0], parts[1]
            qstart, qend   = int(parts[2]), int(parts[3])
            qlen           = int(parts[4])
            sstart, send   = int(parts[5]), int(parts[6])
            strand         = parts[7]
            qseq_aln       = parts[8]
            sseq_aln       = parts[9]

            # For minus-strand hits, BLAST reports sstart > send; SAM uses
            # the lowest coordinate and FLAG 16 for reverse complement
            flag = 0
            read_seq = ref_seqs.get(qseqid, '*')
            if strand == 'minus':
                flag |= 16
                sstart, send = send, sstart
                read_seq = rc(read_seq)

            cigar = pairwise_to_cigar(qseq_aln, sseq_aln, qstart, qend, qlen)
            pos   = sstart   # 1-based leftmost mapping position on reference

            fh.write(
                f'{qseqid}\t{flag}\t{sseqid}\t{pos}\t60\t{cigar}\t*\t0\t0\t'
                f'{read_seq}\t*\n'
            )
            written += 1

print(f"Wrote {written} SAM records to {out_sam}")
PYEOF

    # Convert SAM → sorted BAM
    samtools sort -o "$REFS_BAM" "$REFS_SAM"
    samtools index "$REFS_BAM"
    rm -f "$REFS_SAM"
    echo "[$(date)] Refs BAM written: $REFS_BAM" | tee -a "$LOG"
fi

# Alignment stats + which refs aligned
echo "[$(date)] Reference alignment stats:" | tee -a "$LOG"
samtools flagstat "$REFS_BAM" | tee -a "$LOG"
{
    echo "=== Reference Sequence Alignment Stats ==="
    samtools flagstat "$REFS_BAM"
    echo ""
    echo "Aligned reference accessions and positions on ${CONTIG_ID}:"
    samtools view -F 4 "$REFS_BAM" \
        | awk '{printf "  %-12s  start=%-8d  mapq=%s\n", $1, $4, $5}' \
        | sort -k2 -t= -n
    echo ""
} >> "$REPORT"

echo "[$(date)] Reference positions:" | tee -a "$LOG"
samtools view -F 4 "$REFS_BAM" \
    | awk '{printf "  %-12s  start=%-8d\n", $1, $4}' \
    | sort -k2 -t= -n | tee -a "$LOG"

# ---------------------------------------------------------------------------
# STEP 4: Coverage/pileup at the gap region
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 4: Coverage at gap region (${CONTIG_ID}:${GAP_START}-${GAP_END})"
echo "================================================================="

REGION_START=$(( GAP_START - CONTEXT ))
REGION_END=$(( GAP_END + CONTEXT ))
REGION="${CONTIG_ID}:${REGION_START}-${REGION_END}"

{
    echo "=== Coverage at Gap Region (${CONTIG_ID}:${GAP_START}-${GAP_END}, ±${CONTEXT} bp context) ==="
    echo ""
    printf "%-8s  %-10s  %-10s\n" "POS" "ONT_COV" "REF_COV"
    printf "%-8s  %-10s  %-10s\n" "-------" "-------" "-------"
} | tee -a "$REPORT"

python3 - <<PYEOF 2>&1 | tee -a "$REPORT"
import subprocess

region = "${REGION}"
ont_bam = "${ONT_BAM}"
ref_bam = "${REFS_BAM}"
gap_start = ${GAP_START}
gap_end   = ${GAP_END}
context   = ${CONTEXT}

def get_depth(bam, region):
    """Return dict of pos -> depth for the region."""
    result = subprocess.run(
        ["samtools", "depth", "-a", "-r", region, bam],
        capture_output=True, text=True, check=False
    )
    depth = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            depth[int(parts[1])] = int(parts[2])
    return depth

ont_depth = get_depth(ont_bam, region)
ref_depth = get_depth(ref_bam, region)

all_pos = sorted(set(ont_depth) | set(ref_depth))

in_gap_ont = []
in_gap_ref = []
for pos in all_pos:
    od = ont_depth.get(pos, 0)
    rd = ref_depth.get(pos, 0)
    marker = " <-- GAP" if gap_start <= pos <= gap_end else ""
    print(f"{pos:<8}  {od:<10}  {rd:<10}{marker}")
    if gap_start <= pos <= gap_end:
        in_gap_ont.append(od)
        in_gap_ref.append(rd)

print()
print(f"Gap region summary ({gap_start}-{gap_end}):")
if in_gap_ont:
    print(f"  ONT reads:  min={min(in_gap_ont)}  max={max(in_gap_ont)}  mean={sum(in_gap_ont)/len(in_gap_ont):.1f}")
else:
    print("  ONT reads:  no coverage data")
if in_gap_ref:
    print(f"  References: min={min(in_gap_ref)}  max={max(in_gap_ref)}  mean={sum(in_gap_ref)/len(in_gap_ref):.1f}")
    if max(in_gap_ref) > 0:
        print("  -> At least one reference spans the gap region")
    else:
        print("  -> No reference sequences span the gap region")
else:
    print("  References: no coverage data")
PYEOF

echo "[$(date)] STEP 4 done" | tee -a "$LOG"

# ---------------------------------------------------------------------------
# STEP 5: Extract focused COX1 region for IGV visualization
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 5: Extract COX1 region for IGV visualization"
echo "================================================================="

VIZ_REGION="${CONTIG_ID}:2100-3500"
ONT_REGION_BAM="${OUTDIR}/alignments/cox1_region_ont.bam"
REFS_REGION_BAM="${OUTDIR}/alignments/cox1_region_refs.bam"

# Index the assembly FASTA so IGV can use it as the genome reference directly.
# BAM coordinates reference the full contig, so the full assembly must be the
# IGV genome — a region-extracted FASTA has the wrong coordinate space.
if [[ ! -s "${ASSEMBLY}.fai" ]]; then
    samtools faidx "$ASSEMBLY"
fi

# Subset BAMs to COX1 region — regenerate if source BAM is newer
if [[ ! -s "$ONT_REGION_BAM" || "$ONT_BAM" -nt "$ONT_REGION_BAM" ]]; then
    samtools view -b -h "$ONT_BAM" "$VIZ_REGION" \
        | samtools sort -o "$ONT_REGION_BAM"
    samtools index "$ONT_REGION_BAM"
    echo "[$(date)] ONT region BAM: $ONT_REGION_BAM" | tee -a "$LOG"
fi

if [[ ! -s "$REFS_REGION_BAM" || "$REFS_BAM" -nt "$REFS_REGION_BAM" ]]; then
    samtools view -b -h "$REFS_BAM" "$VIZ_REGION" \
        | samtools sort -o "$REFS_REGION_BAM"
    samtools index "$REFS_REGION_BAM"
    echo "[$(date)] Refs region BAM: $REFS_REGION_BAM" | tee -a "$LOG"
fi

echo "[$(date)] STEP 5 done" | tee -a "$LOG"

# ---------------------------------------------------------------------------
# STEP 6: Summary report
# ---------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "STEP 6: Summary"
echo "================================================================="

python3 - <<PYEOF 2>&1 | tee -a "$REPORT"
import subprocess

ont_bam  = "${ONT_BAM}"
ref_bam  = "${REFS_BAM}"
gap_s    = ${GAP_START}
gap_e    = ${GAP_END}
contig   = "${CONTIG_ID}"
region   = f"{contig}:{gap_s}-{gap_e}"

# Count refs that span the gap
spanning = []
result = subprocess.run(
    ["samtools", "view", "-F", "4", ref_bam],
    capture_output=True, text=True, check=False
)
for line in result.stdout.strip().split("\n"):
    if not line:
        continue
    parts = line.split("\t")
    if len(parts) < 10:
        continue
    acc   = parts[0]
    start = int(parts[3])
    end   = start + len(parts[9]) - 1
    if start <= gap_s and end >= gap_e:
        spanning.append((acc, start, end))

print("\n=== Summary ===")
print(f"Gap region: {contig}:{gap_s}-{gap_e} (13 bp, sequence: CCCTCACCCCCTC)")
print(f"References spanning the gap: {len(spanning)}")
for acc, s, e in spanning:
    print(f"  {acc}  ({s}-{e})")

if spanning:
    print()
    print("RECOMMENDATION: Reference sequences cover the gap. Extract the reference")
    print("  sequence at this position to determine the correct bases, then decide")
    print("  whether to correct the assembly or submit with a note.")
else:
    print()
    print("RECOMMENDATION: No reference sequences span the gap. The 13 bp region")
    print("  may be outside the amplified barcode region for these references.")
    print("  Examine ONT read coverage to assess whether the gap is a real deletion.")

print()
print("IGV files for inspection:")
print(f"  Genome (full assembly): assemblies/mtDNA_02/assembly.fasta")
print(f"  ONT reads track:        alignments/cox1_region_ont.bam")
print(f"  Reference seqs track:   alignments/cox1_region_refs.bam")
print(f"  Navigate to:            {contig}:{gap_s - 100}-{gap_e + 100}")
PYEOF

echo "" | tee -a "$REPORT"
echo "Full report: $REPORT" | tee -a "$LOG"
echo ""
echo "[$(date)] align_cox1_locus.sh complete." | tee -a "$LOG"
echo ""
echo "To inspect in IGV:"
echo "  1. Genome (full assembly): ${ASSEMBLY}"
echo "  2. Track 1 (ONT reads):    ${OUTDIR}/alignments/cox1_region_ont.bam"
echo "  3. Track 2 (references):   ${OUTDIR}/alignments/cox1_region_refs.bam"
echo ""
echo "  Navigate to: ${CONTIG_ID}:$(( GAP_START - 100 ))-$(( GAP_END + 100 ))"
