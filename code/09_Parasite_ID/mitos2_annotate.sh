#!/usr/bin/env bash
# =============================================================================
# MITOS2 Local Annotation Script
# Dioctophyme renale Mitochondrial Minichromosomes
# =============================================================================
# Annotates each assembled minichromosome with MITOS2, matching the Galaxy
# server settings used in Macchiaroli et al. 2025:
#   Genetic code : Invertebrate Mitochondrial (5)
#   Reference    : RefSeq89 Metazoa
#
# Dependencies:
#   mitos (bioconda)  — conda install -c bioconda mitos
#   Reference data    — downloaded automatically on first run
#     DOI: 10.5281/zenodo.3685310
#
# Usage:
#   bash mitos2_annotate.sh
#
# Outputs:
#   mt_annotations/mtDNA_XX/result.gff   — gene coordinates
#   mt_annotations/mtDNA_XX/result.bed   — BED format
#   mt_annotations/mtDNA_XX/result.mitos — annotated feature sequences
#   mt_annotations/mtDNA_XX/result.fas   — full annotated chromosome
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration — edit paths here if needed
# =============================================================================

GENETIC_CODE=5                              # 5 = Invertebrate Mitochondrial
REFDATA_DIR="mitos2_refdata"               # Parent dir for reference data
REFSEQ_DIR="${REFDATA_DIR}/refseq89m"      # RefSeq89 Metazoa subdirectory
REFDATA_URL="https://zenodo.org/records/3685310/files/refseq89m.tar.bz2"
ASSEMBLY_DIR="mt_assemblies"
ANNOTATION_DIR="mt_annotations"

# =============================================================================
# STEP 0: Preflight checks
# =============================================================================

echo "=== Checking MITOS2 installation ==="
if ! command -v runmitos.py &>/dev/null; then
  cat <<'EOF'
ERROR: runmitos.py not found in PATH.

Install into a new environment:
  conda create -n mitos2 -c bioconda -c conda-forge mitos
  conda activate mitos2

Or install into the current environment:
  conda install -c bioconda mitos

Then re-run this script with the mitos2 environment active.
EOF
  exit 1
fi
echo "Found: $(which runmitos.py)"

echo ""
echo "=== Checking RefSeq89 Metazoa reference data ==="
if [ ! -d "$REFSEQ_DIR" ]; then
  echo "Reference data not found at: ${REFSEQ_DIR}"
  echo "Downloading from Zenodo (DOI: 10.5281/zenodo.3685310) — ~20MB..."
  mkdir -p "$REFDATA_DIR"
  curl -L --progress-bar \
    "$REFDATA_URL" \
    -o "${REFDATA_DIR}/refseq89m.tar.bz2"
  echo "Extracting..."
  tar -xjf "${REFDATA_DIR}/refseq89m.tar.bz2" -C "$REFDATA_DIR/"
  rm "${REFDATA_DIR}/refseq89m.tar.bz2"
  echo "Reference data ready: ${REFSEQ_DIR}"
else
  echo "Found: ${REFSEQ_DIR}"
fi

mkdir -p "$ANNOTATION_DIR"

# =============================================================================
# STEP 1: Annotate each minichromosome
# =============================================================================
# MITOS2 flags (matching Galaxy server settings):
#   -i  input FASTA
#   -c  genetic code (5 = invertebrate mitochondrial)
#   -o  output directory
#   -R  base directory containing all reference data (--refdir)
#   -r  subdirectory name within -R for this reference version (--refseqver)
#       MITOS2 joins these internally: path = REFDIR/REFSEQVER
#       so -R mitos2_refdata -r refseq89m → mitos2_refdata/refseq89m
#   --noplots  skip R-based circular genome plots (avoids R dependency)

echo ""
echo "=== Running MITOS2 annotation ==="

for gene in 01 02 03 04 05 06 07 08 09 10 11; do
  asm="${ASSEMBLY_DIR}/mtDNA_${gene}/assembly.fasta"
  outdir="${ANNOTATION_DIR}/mtDNA_${gene}"

  if [ ! -f "$asm" ] || [ ! -s "$asm" ]; then
    echo ""
    echo "mtDNA_${gene}: no assembly — skipping"
    continue
  fi

  n_seqs=$(grep -c "^>" "$asm")
  asm_len=$(awk '!/^>/{sum+=length($0)} END{print sum}' "$asm")
  echo ""
  echo "mtDNA_${gene}: ${n_seqs} contig(s), ${asm_len}bp"

  # MITOS2 requires the output directory to not already exist for this run.
  # We remove and recreate so re-runs don't accumulate stale results.
  rm -rf "$outdir"
  mkdir -p "$outdir"

  runmitos.py \
    -i "$asm" \
    -c "$GENETIC_CODE" \
    -o "$outdir" \
    -R "$REFDATA_DIR" \
    -r "refseq89m" \
    --noplots \
    2>&1 | grep -vE "^(DEBUG|WARNING: No handlers)" || true

  # Report results
  gff="${outdir}/result.gff"
  if [ -f "$gff" ] && [ -s "$gff" ]; then
    n_pcg=$(grep -c "gene_id" "$gff" 2>/dev/null \
      | head -1 || echo 0)
    genes=$(grep "Name=" "$gff" 2>/dev/null \
      | sed 's/.*Name=\([^;]*\).*/\1/' | sort -u | tr '\n' ' ' || true)
    echo "  Features: ${genes}"
  else
    echo "  WARNING: no result.gff produced — check ${outdir}/ for error output"
  fi
done

# =============================================================================
# STEP 2: Summary table
# =============================================================================

echo ""
echo "=== Annotation summary ==="
echo ""
printf "%-12s %-10s %-6s  %s\n" "Chromosome" "Length(bp)" "Status" "Annotated genes"
printf "%-12s %-10s %-6s  %s\n" "----------" "----------" "------" "---------------"

for gene in 01 02 03 04 05 06 07 08 09 10 11; do
  asm="${ASSEMBLY_DIR}/mtDNA_${gene}/assembly.fasta"
  gff="${ANNOTATION_DIR}/mtDNA_${gene}/result.gff"

  if [ ! -f "$asm" ] || [ ! -s "$asm" ]; then
    printf "%-12s %-10s %-6s  %s\n" "mtDNA_${gene}" "N/A" "SKIP" "no assembly"
    continue
  fi

  asm_len=$(awk '!/^>/{sum+=length($0)} END{print sum}' "$asm")

  if [ ! -f "$gff" ] || [ ! -s "$gff" ]; then
    printf "%-12s %-10s %-6s  %s\n" "mtDNA_${gene}" "$asm_len" "FAIL" "no annotation"
    continue
  fi

  genes=$(grep "Name=" "$gff" 2>/dev/null \
    | sed 's/.*Name=\([^;]*\).*/\1/' | sort -u | tr '\n' ' ' || true)
  printf "%-12s %-10s %-6s  %s\n" "mtDNA_${gene}" "$asm_len" "OK" "$genes"
done

echo ""
echo "=== Complete ==="
echo "Full results in: ${ANNOTATION_DIR}/"
echo ""
echo "Key files per chromosome:"
echo "  result.gff   — gene coordinates for genome browsers (IGV, etc.)"
echo "  result.bed   — BED format, easy to intersect/filter"
echo "  result.mitos — MITOS feature sequences (FASTA)"
echo "  result.fas   — full annotated chromosome sequence"
echo ""
echo "Next: manual curation in IGV"
echo "  minimap2 --map-ont assembly.fasta reads.fa | samtools sort > aln.bam"
echo "  samtools index aln.bam"
echo "  # Open IGV: load assembly as genome, load aln.bam and result.gff"
echo "  # Translation table: Invertebrate Mitochondrial"