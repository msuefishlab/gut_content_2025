#!/usr/bin/env bash
# =============================================================================
# identify_worm_sequences.sh
# Pipeline to identify taxonomically informative sequences in ONT reads
# using Barrnap (rRNA), minimap2 vs BOLD (COI), and BLASTn vs NCBI nr.
#
# Usage:
#   bash identify_worm_sequences.sh -b BOLD_FASTA [options] READS1.fa READS2.fa ...
#
# Options:
#   -b  PATH   BOLD FASTA database (required)
#   -o  PATH   Output directory [default: worm_id_results]
#   -k  STR    Barrnap kingdom: euk|bac|mito [default: euk]
#   -n  PATH   Local BLAST nr database path (e.g. /db/nr/nr)
#              If omitted, BLASTn runs remotely against NCBI nr (-remote)
#   -e  FLOAT  BLASTn e-value cutoff [default: 1e-10]
#   -t  INT    BLASTn max target sequences [default: 10]
#   -T  INT    Threads for minimap2 and local BLASTn [default: 4]
#   -h         Show this help
#
# Dependencies: barrnap, samtools, minimap2, blastn, python3 (Biopython)
#
# Example (remote BLAST):
#   bash identify_worm_sequences.sh \
#       -b BOLD_Public.27-Mar-2026.fasta \
#       -o results \
#       bc1_reads.fa bc2_reads.fa
#
# Example (local nr):
#   bash identify_worm_sequences.sh \
#       -b BOLD_Public.27-Mar-2026.fasta \
#       -n /databases/nr/nr \
#       -T 16 \
#       -o results \
#       bc1_reads.fa bc2_reads.fa
# =============================================================================

set -euo pipefail

# --- Defaults ----------------------------------------------------------------
OUTDIR="worm_id_results"
BOLD_FA=""
KINGDOM="euk"
NR_DB=""          # empty = use -remote
EVALUE="1e-10"
MAX_SEQS=10
THREADS=4

# --- Argument parsing --------------------------------------------------------
usage() {
    sed -n '2,28p' "$0" | sed 's/^# //'
    exit 1
}

while getopts "b:o:k:n:e:t:T:h" opt; do
    case $opt in
        b) BOLD_FA="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        k) KINGDOM="$OPTARG" ;;
        n) NR_DB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) MAX_SEQS="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

READS=("$@")

if [[ -z "$BOLD_FA" || ${#READS[@]} -eq 0 ]]; then
    echo "ERROR: Must provide -b BOLD_FASTA and at least one reads file."
    usage
fi

# --- Setup -------------------------------------------------------------------
mkdir -p "$OUTDIR"/{barrnap,bold,blast_queries,blast_results}
LOG="$OUTDIR/pipeline.log"
echo "[$(date)] Pipeline started" | tee "$LOG"

if [[ -n "$NR_DB" ]]; then
    echo "[$(date)] BLASTn mode: local nr database: $NR_DB" | tee -a "$LOG"
else
    echo "[$(date)] BLASTn mode: remote NCBI nr (--remote)" | tee -a "$LOG"
    echo "[$(date)] NOTE: remote BLAST can be slow; consider -n for local nr if available." | tee -a "$LOG"
fi

# =============================================================================
# STEP 1: Filter BOLD to Nematoda
# =============================================================================
BOLD_NEMA="$OUTDIR/bold/bold_nematoda.fa"
BOLD_MMI="$OUTDIR/bold/bold_nematoda.mmi"

if [[ ! -f "$BOLD_NEMA" ]]; then
    echo "[$(date)] Filtering BOLD to Nematoda..." | tee -a "$LOG"
    python3 - <<PYEOF
from Bio import SeqIO
count = 0
with open("$BOLD_NEMA", "w") as out:
    for rec in SeqIO.parse("$BOLD_FA", "fasta"):
        if "Nematoda" in rec.description:
            SeqIO.write(rec, out, "fasta")
            count += 1
print(f"  Wrote {count} Nematoda sequences to $BOLD_NEMA")
PYEOF
else
    echo "[$(date)] BOLD Nematoda subset already exists, skipping filter." | tee -a "$LOG"
fi

NEMA_COUNT=$(grep -c ">" "$BOLD_NEMA")
echo "[$(date)] Nematoda BOLD sequences: $NEMA_COUNT" | tee -a "$LOG"

# Build minimap2 index
if [[ ! -f "$BOLD_MMI" ]]; then
    echo "[$(date)] Building minimap2 index for BOLD Nematoda..." | tee -a "$LOG"
    minimap2 -d "$BOLD_MMI" "$BOLD_NEMA" 2>>"$LOG"
fi

# =============================================================================
# STEP 2: Process each sample
# =============================================================================
for READS_FA in "${READS[@]}"; do

    SAMPLE=$(basename "$READS_FA" | sed 's/\.[^.]*$//')
    echo "" | tee -a "$LOG"
    echo "[$(date)] ========== Processing sample: $SAMPLE ==========" | tee -a "$LOG"

    SDIR_BARRNAP="$OUTDIR/barrnap/$SAMPLE"
    SDIR_BOLD="$OUTDIR/bold/$SAMPLE"
    mkdir -p "$SDIR_BARRNAP" "$SDIR_BOLD"

    # -------------------------------------------------------------------------
    # STEP 2a: Index reads for extraction
    # -------------------------------------------------------------------------
    if [[ ! -f "${READS_FA}.fai" ]]; then
        echo "[$(date)]   Indexing $READS_FA..." | tee -a "$LOG"
        samtools faidx "$READS_FA"
    fi

    # -------------------------------------------------------------------------
    # STEP 2b: Barrnap rRNA detection
    # -------------------------------------------------------------------------
    GFF="$SDIR_BARRNAP/${SAMPLE}_euk.gff"
    echo "[$(date)]   Running Barrnap (--kingdom $KINGDOM)..." | tee -a "$LOG"
    barrnap --kingdom "$KINGDOM" "$READS_FA" > "$GFF" 2>>"$LOG"

    NFEAT=$(grep -v "^#" "$GFF" | wc -l)
    echo "[$(date)]   Barrnap found $NFEAT rRNA features." | tee -a "$LOG"

    # -------------------------------------------------------------------------
    # STEP 2c: Extract annotated rRNA regions + putative ITS1
    # -------------------------------------------------------------------------
    echo "[$(date)]   Extracting rRNA sequences..." | tee -a "$LOG"

    python3 - <<PYEOF
import subprocess, os
from collections import defaultdict

gff_file = "$GFF"
reads_fa = "$READS_FA"
out_dir  = "$OUTDIR/blast_queries"
sample   = "$SAMPLE"
log_file = "$LOG"

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
            "seq_id": seq_id, "type": feat_type, "name": name,
            "start": int(start), "end": int(end),
            "strand": strand, "score": score
        })

by_read = defaultdict(list)
for f in features:
    by_read[f["seq_id"]].append(f)

extracted = []

for seq_id, feats in by_read.items():
    # Extract each individual rRNA feature (skip 5S)
    for f in feats:
        if "5S_rRNA" in f["name"]:
            continue
        region = f"{seq_id}:{f['start']}-{f['end']}"
        tag = f"{sample}_{f['name']}_{seq_id[:8]}"
        outfa = os.path.join(out_dir, f"{tag}.fa")
        result = subprocess.run(["samtools", "faidx", reads_fa, region],
                                capture_output=True, text=True)
        fasta = result.stdout.replace(region, tag + f"|{f['name']}|{f['score']}", 1)
        with open(outfa, "w") as fh:
            fh.write(fasta)
        extracted.append(outfa)
        msg = f"  Extracted {f['name']} from {seq_id[:8]} -> {outfa}"
        print(msg)
        open(log_file, "a").write(msg + "\n")

    # If same read has both 18S and 5.8S, extract intervening ITS1 region
    names = [f["name"] for f in feats]
    if "18S_rRNA" in names and "5_8S_rRNA" in names:
        f18 = next(f for f in feats if f["name"] == "18S_rRNA")
        f58 = next(f for f in feats if f["name"] == "5_8S_rRNA")
        its_start = f18["end"] + 1
        its_end   = f58["start"] - 1
        if its_end > its_start:
            region = f"{seq_id}:{its_start}-{its_end}"
            tag = f"{sample}_ITS1_putative_{seq_id[:8]}"
            outfa = os.path.join(out_dir, f"{tag}.fa")
            result = subprocess.run(["samtools", "faidx", reads_fa, region],
                                    capture_output=True, text=True)
            fasta = result.stdout.replace(region, tag + "|ITS1_putative|extracted", 1)
            with open(outfa, "w") as fh:
                fh.write(fasta)
            extracted.append(outfa)
            msg = f"  Extracted putative ITS1 from {seq_id[:8]} ({its_start}-{its_end}) -> {outfa}"
            print(msg)
            open(log_file, "a").write(msg + "\n")

print(f"  Total sequences extracted: {len(extracted)}")
PYEOF

    # -------------------------------------------------------------------------
    # STEP 2d: minimap2 vs BOLD Nematoda
    # -------------------------------------------------------------------------
    BAM="$SDIR_BOLD/${SAMPLE}_vs_bold_nematoda.bam"
    HITS_TXT="$SDIR_BOLD/${SAMPLE}_bold_hits.txt"

    echo "[$(date)]   Aligning $SAMPLE reads vs. BOLD Nematoda..." | tee -a "$LOG"
    minimap2 -ax map-ont \
        --secondary=no \
        -f 0.0001 \
        -t "$THREADS" \
        "$BOLD_MMI" \
        "$READS_FA" \
        2>>"$LOG" \
    | samtools view -F 4 -b \
    | samtools sort -o "$BAM"
    samtools index "$BAM"

    echo "[$(date)]   BOLD alignment summary for $SAMPLE:" | tee -a "$LOG"
    samtools flagstat "$BAM" | tee -a "$LOG"

    echo "# read_id  reference  mapq  flag" > "$HITS_TXT"
    samtools view "$BAM" | awk '{print $1, $3, $5, $2}' >> "$HITS_TXT"

    echo "[$(date)]   Top BOLD hits (by mapq):" | tee -a "$LOG"
    grep -v "^#" "$HITS_TXT" | sort -k3 -rn | head -10 | tee -a "$LOG"

done

# =============================================================================
# STEP 3: BLASTn extracted rRNA sequences against nr
# =============================================================================
echo "" | tee -a "$LOG"
echo "[$(date)] ========== BLASTn rRNA sequences vs. NCBI nr ==========" | tee -a "$LOG"

# Tabular output fields:
# qseqid sseqid pident length qcovs evalue bitscore stitle sscinames
BLAST_FMT="6 qseqid sseqid pident length qcovs evalue bitscore stitle sscinames"

shopt -s nullglob
QUERY_FILES=("$OUTDIR/blast_queries/"*.fa)

if [[ ${#QUERY_FILES[@]} -eq 0 ]]; then
    echo "[$(date)] WARNING: No query sequences found — skipping BLASTn." | tee -a "$LOG"
else
    echo "[$(date)] Found ${#QUERY_FILES[@]} query sequence(s) to BLAST." | tee -a "$LOG"

    for QUERY_FA in "${QUERY_FILES[@]}"; do
        QNAME=$(basename "$QUERY_FA" .fa)
        BLAST_TSV="$OUTDIR/blast_results/${QNAME}.blast.tsv"
        BLAST_TXT="$OUTDIR/blast_results/${QNAME}.blast.txt"

        echo "[$(date)]   BLASTn: $QNAME" | tee -a "$LOG"

        if [[ -n "$NR_DB" ]]; then
            # --- Local nr --------------------------------------------------------
            blastn \
                -query    "$QUERY_FA" \
                -db       "$NR_DB" \
                -out      "$BLAST_TSV" \
                -outfmt   "$BLAST_FMT" \
                -evalue   "$EVALUE" \
                -max_target_seqs "$MAX_SEQS" \
                -num_threads "$THREADS" \
                2>>"$LOG"

            # Human-readable pairwise alignment (top 5 only)
            blastn \
                -query    "$QUERY_FA" \
                -db       "$NR_DB" \
                -out      "$BLAST_TXT" \
                -outfmt   0 \
                -evalue   "$EVALUE" \
                -max_target_seqs 5 \
                -num_threads "$THREADS" \
                2>>"$LOG"
        else
            # --- Remote NCBI BLAST -----------------------------------------------
            # NOTE: -remote does not support -num_threads
            blastn \
                -query    "$QUERY_FA" \
                -db       nr \
                -out      "$BLAST_TSV" \
                -outfmt   "$BLAST_FMT" \
                -evalue   "$EVALUE" \
                -max_target_seqs "$MAX_SEQS" \
                -remote \
                2>>"$LOG"

            blastn \
                -query    "$QUERY_FA" \
                -db       nr \
                -out      "$BLAST_TXT" \
                -outfmt   0 \
                -evalue   "$EVALUE" \
                -max_target_seqs 5 \
                -remote \
                2>>"$LOG"
        fi

        # Report top hit to log
        if [[ -s "$BLAST_TSV" ]]; then
            echo "[$(date)]     Top hit:" | tee -a "$LOG"
            head -1 "$BLAST_TSV" | awk -F'\t' \
                '{printf "     %s -> %s | pident=%.1f%% qcov=%s%% evalue=%s\n", \
                  $1, $9, $3, $5, $6}' \
                | tee -a "$LOG"
        else
            echo "[$(date)]     No hits above evalue $EVALUE." | tee -a "$LOG"
        fi
    done
fi

# =============================================================================
# STEP 4: Summary
# =============================================================================
echo "" | tee -a "$LOG"
echo "[$(date)] ========== SUMMARY ==========" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "--- BOLD COI hits ---" | tee -a "$LOG"
for f in "$OUTDIR"/bold/*/*.txt; do
    [[ -f "$f" ]] || continue
    echo "  $(basename "$(dirname "$f")"):" | tee -a "$LOG"
    grep -v "^#" "$f" | awk '{print "    ", $2, "mapq="$3}' | sort -u | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- BLASTn top hits (rRNA/ITS) ---" | tee -a "$LOG"
for f in "$OUTDIR/blast_results/"*.blast.tsv; do
    [[ -s "$f" ]] || continue
    echo "  $(basename "$f" .blast.tsv):" | tee -a "$LOG"
    head -3 "$f" | awk -F'\t' \
        '{printf "    %s | pident=%.1f%% qcov=%s%% evalue=%s\n", $9, $3, $5, $6}' \
        | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "[$(date)] Done." | tee -a "$LOG"
echo "Full log:           $LOG"
echo "rRNA BLAST queries: $OUTDIR/blast_queries/"
echo "BLASTn results:     $OUTDIR/blast_results/  (.tsv = tabular, .txt = pairwise)"
echo "BOLD COI results:   $OUTDIR/bold/"
