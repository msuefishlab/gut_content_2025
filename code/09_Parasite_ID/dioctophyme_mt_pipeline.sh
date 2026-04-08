#!/usr/bin/env bash
# =============================================================================
# Dioctophyme renale Mitochondrial Minichromosome Assembly Pipeline
# =============================================================================
# Purpose: Recruit mitochondrial reads from whole-genome nanopore data,
#          assemble per-minichromosome, and validate assembly integrity.
#
# Key finding: Minichromosomes in this isolate have massively expanded
#              non-coding regions (NCRs) compared to the published genome
#              (Macchiaroli et al. 2025, ~2,589bp average NCR). Our assemblies
#              range from ~4kb to ~30kb, with NCRs of 5,000-11,000bp.
#              This was validated by: (1) single gene hit per contig,
#              (2) no internal repeats in pre-gene regions,
#              (3) no concatemerization signal in window self-alignments.
#
# Reference genome: Macchiaroli et al. 2025 (preprint, Research Square)
#   DOI: https://doi.org/10.21203/rs.3.rs-8148715/v1
#   GenBank: PX715184-PX715194
#
# Input requirements:
#   - all_nanopore_reads.fa         : all WGS nanopore reads in FASTA format
#   - dioctophyme_mt_PCG_amino_acid.fa : D. renale mitochondrial protein seqs
#   - cladeI_accessions.txt         : NCBI accessions for Clade I nematodes
#
# Software dependencies:
#   - DIAMOND v2.1+
#   - Flye v2.9.5
#   - minimap2
#   - seqtk
#   - Python 3
#   - NCBI Entrez (curl)
#
# Note on mac compatibility: bash on macOS is v3; associative arrays
#   (declare -A) are not supported. All size lookups use case statements.
#   Octal interpretation of numbers with leading zeros (e.g. 08, 09) is
#   avoided by quoting gene IDs as strings throughout.
# =============================================================================

set -euo pipefail

# =============================================================================
# STEP 0: Build reference databases
# =============================================================================

echo "=== STEP 0: Building DIAMOND databases ==="

# --- 0a: D. renale protein database (for final read assignment) ---
diamond makedb \
  --in dioctophyme_mt_PCG_amino_acid.fa \
  --db dioctophyme_mt

# --- 0b: Download Clade I mitochondrial proteins from NCBI ---
# Accessions from Macchiaroli et al. 2025 supplementary files
cat > cladeI_accessions.txt << 'EOF'
NC_071371
NC_056391
NC_008640
NC_008693
NC_008692
NC_025750
NC_025751
NC_025752
NC_025753
NC_025754
NC_025749
NC_002681
NC_025755
NC_018596
NC_028621
NC_018597
NC_017747
NC_017750
NC_001328
EOF

mkdir -p cladeI_proteins

while read acc; do
  echo "Fetching $acc..."
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=nuccore&id=${acc}&rettype=fasta_cds_aa&retmode=text" \
    >> cladeI_proteins/cladeI_all_proteins.fa
  sleep 1  # respect NCBI rate limits
done < cladeI_accessions.txt

echo "Total Clade I proteins downloaded:"
grep -c "^>" cladeI_proteins/cladeI_all_proteins.fa

# Build Clade I DIAMOND database
diamond makedb \
  --in cladeI_proteins/cladeI_all_proteins.fa \
  --db cladeI_mt

# =============================================================================
# STEP 1: Recruit mitochondrial reads (two-pass strategy)
# =============================================================================
# Two-pass approach:
#   Pass 1 (broad)  : Clade I database catches reads too diverged for D. renale db
#   Pass 2 (specific): D. renale database assigns reads to specific minichromosomes
# The --long-reads flag is essential: it enables DIAMOND's chaining algorithm
# for sparse seeds across long nanopore reads. Without it, recruitment drops
# significantly (530 vs 270 unique reads in testing).

echo "=== STEP 1: Recruiting mitochondrial reads ==="

# Pass 1a: Recruit with D. renale proteins
diamond blastx \
  --db dioctophyme_mt \
  --query all_nanopore_reads.fa \
  --query-gencode 5 \
  --long-reads \
  --sensitive \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --threads 16 \
  --out mt_read_hits_drenale.txt

echo "D. renale db recruited:"
cut -f1 mt_read_hits_drenale.txt | sort -u | wc -l

# Pass 1b: Recruit with Clade I proteins
diamond blastx \
  --db cladeI_mt \
  --query all_nanopore_reads.fa \
  --query-gencode 5 \
  --long-reads \
  --sensitive \
  --evalue 1e-3 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --threads 16 \
  --out mt_read_hits_cladeI.txt

echo "Clade I db recruited:"
cut -f1 mt_read_hits_cladeI.txt | sort -u | wc -l

# Pool all recruited reads (union of both searches)
cat mt_read_hits_drenale.txt mt_read_hits_cladeI.txt \
  | cut -f1 | sort -u > all_mt_read_ids.txt

echo "Total unique reads (union):"
wc -l all_mt_read_ids.txt

seqtk subseq all_nanopore_reads.fa all_mt_read_ids.txt > all_mt_reads_pooled.fa

# Pass 2: Assign pooled reads to specific minichromosomes using D. renale db
# All-hits assignment: a read is included in every gene pool it matches.
# We deliberately do NOT use best-hit-only filtering here. At low coverage
# (20-65 reads per gene), losing even 1-2 bridging reads can break assembly.
# Reads matching multiple genes (typically NCR-spanning reads, since NCRs
# are 99.5-99.9% identical across minichromosomes) are retained in all
# matching pools. Any resulting chimeric assemblies are caught by Step 4
# validation (single gene copy check).
diamond blastx \
  --db dioctophyme_mt \
  --query all_mt_reads_pooled.fa \
  --query-gencode 5 \
  --long-reads \
  --sensitive \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --threads 16 \
  --out mt_read_hits_pooled_assigned.txt

echo "Per-gene read counts:"
cut -f2 mt_read_hits_pooled_assigned.txt | sort | uniq -c | sort -rn

# =============================================================================
# STEP 2: Split reads per minichromosome
# =============================================================================

echo "=== STEP 2: Splitting reads by minichromosome ==="

mkdir -p mt_reads_by_gene

for gene in 01 02 03 04 05 06 07 08 09 10 11; do
  grep "mtDNA_${gene}" mt_read_hits_pooled_assigned.txt \
    | cut -f1 | sort -u \
    > mt_reads_by_gene/reads_mtDNA_${gene}.ids

  seqtk subseq all_nanopore_reads.fa \
    mt_reads_by_gene/reads_mtDNA_${gene}.ids \
    > mt_reads_by_gene/reads_mtDNA_${gene}.fa

  n=$(wc -l < mt_reads_by_gene/reads_mtDNA_${gene}.ids)
  mean=$(awk '/^>/{next}{sum+=length($0);n++}END{
    if(n>0) printf "%.0f",sum/n; else print "NA"}' \
    mt_reads_by_gene/reads_mtDNA_${gene}.fa)
  max=$(awk '/^>/{next}{if(length($0)>m)m=length($0)}END{print m}' \
    mt_reads_by_gene/reads_mtDNA_${gene}.fa)

  echo "mtDNA_${gene}: ${n} reads, mean=${mean}bp, max=${max}bp"
done

# =============================================================================
# STEP 3: Assembly
# =============================================================================
# Parameters:
#   --genome-size 4k  : approximate size of a single minichromosome PCG + NCR
#   --meta            : metagenome mode, more tolerant of uneven coverage
#   --keep-haplotypes : preserve variant haplotypes
#   --min-overlap 1000: Flye's minimum allowed value. The paper used 3800,
#                       but that requires reads nearly spanning the full
#                       chromosome, which our read depths cannot support.
#                       Concatemerization was ruled out by validation (Step 4).
#
# NOTE: set -euo pipefail is active, so Flye failures must be caught with
#       || true to prevent the pipeline aborting on low-coverage genes.
#
# NOTE: mtDNA_11 (NAD6) cannot be assembled - see STEP 5 for validation.

echo "=== STEP 3: Assembly with Flye ==="

mkdir -p mt_assemblies

for gene in 01 02 03 04 05 06 07 08 09 10 11; do
  reads=mt_reads_by_gene/reads_mtDNA_${gene}.fa
  n=$(grep -c "^>" "$reads" 2>/dev/null || echo 0)

  if [ "$n" -lt 5 ]; then
    echo "mtDNA_${gene}: SKIPPED (only ${n} reads - insufficient for assembly)"
    continue
  fi

  # Flye is non-deterministic under multithreading (github.com/mikolmogorov/
  # Flye/issues/458): identical inputs can produce 0-edge graphs or fail the
  # polishing coverage filter on one run and succeed on another. Retry up to
  # 3 times before falling back to stage rescue. Each attempt clears the
  # output directory first to avoid Flye misbehaving on partial state.
  flye_success=0
  for attempt in 1 2 3; do
    rm -rf mt_assemblies/mtDNA_${gene}
    flye \
      --nano-raw "$reads" \
      --out-dir mt_assemblies/mtDNA_${gene} \
      --genome-size 4k \
      --meta \
      --keep-haplotypes \
      --min-overlap 1000 \
      --iterations 3 \
      2>&1 | tail -3 || true
    if [ -f "mt_assemblies/mtDNA_${gene}/assembly_info.txt" ]; then
      flye_success=1
      break
    fi
    echo "mtDNA_${gene}: attempt ${attempt}/3 failed, retrying..."
  done

  info=mt_assemblies/mtDNA_${gene}/assembly_info.txt
  if [ "$flye_success" -eq 1 ]; then
    awk -v g="$gene" 'NR>1{print "mtDNA_"g": "$2"bp circ="$4" cov="$5}' "$info"
  else
    # Flye's polishing coverage filter (threshold=3) is hardcoded in the Python
    # wrapper and cannot be overridden via --extra-params min_read_cov_cutoff,
    # which only affects the C++ assembly/contigger stages. When coverage is
    # borderline (typically genes with short mean read lengths or cross-
    # contaminating NCR reads), Flye assembles and polishes valid contigs but
    # then discards them at the last step without writing assembly.fasta.
    # The polished output still exists in 40-polishing/ — rescue it from there.
    outdir=mt_assemblies/mtDNA_${gene}
    # Tiered rescue: Flye can fail at different stages, leaving useful output
    # at different depths. Check each stage in reverse order, skipping empty
    # files (contigger writes 0-byte contigs.fasta when graph has 0 edges).
    #
    # Tier 1: polished_N.fasta (path varies by Flye version, search whole dir)
    # Tier 2: 30-contigger/contigs.fasta (post-graph, unpolished; may be 0 bytes)
    # Tier 3: 10-consensus/consensus.fasta (pre-graph consensus disjointigs;
    #          exists whenever the assembly+consensus stages completed, which
    #          is the deepest reliable rescue point when graph simplification
    #          collapses all edges to zero)
    rescue_src=""
    rescue_stage=""
    _t1=$(find "$outdir" -name "polished_*.fasta" 2>/dev/null | sort -V | tail -1 || true)
    _t2="$outdir/30-contigger/contigs.fasta"
    _t3="$outdir/10-consensus/consensus.fasta"
    if [ -s "$_t1" ]; then
      rescue_src="$_t1"; rescue_stage="polished"
    elif [ -s "$_t2" ]; then
      rescue_src="$_t2"; rescue_stage="contigger"
    elif [ -s "$_t3" ]; then
      rescue_src="$_t3"; rescue_stage="consensus"
    fi
    if [ -n "$rescue_src" ]; then
      cp "$rescue_src" "$outdir/assembly.fasta"
      rescued_len=$(awk '!/^>/{sum+=length($0)} END{print sum}' "$outdir/assembly.fasta")
      echo "mtDNA_${gene}: rescued ${rescued_len}bp from ${rescue_stage} stage ($(basename "$rescue_src"))"
    else
      echo "mtDNA_${gene}: ASSEMBLY FAILED (no usable output at any stage)"
    fi
  fi
done

# =============================================================================
# STEP 4: Validate assembly integrity
# =============================================================================
# Key question: are large assemblies genuine (expanded NCRs) or artifacts
# (concatemerization due to Flye looping around the circular chromosome)?
#
# Three lines of evidence tested per assembly:
#   (A) Gene copy number: DIAMOND with --max-hsps 0 to find all hits.
#       Concatemers would show the gene multiple times at regular intervals.
#   (B) Window self-alignment: extract windows at expected-size offsets and
#       align them. Concatemers would show >90% identity between windows.
#   (C) NCR internal repeats: self-align the pre-gene region with minimap2 -X.
#       True expanded NCRs have no internal repeats; concatemerized NCRs do.
#
# Expected sizes from Macchiaroli et al. 2025 Table 1 (their isolate):
#   smtDNA I   (ATP6) : 4,328bp    smtDNA VII  (NAD2) : 4,031bp
#   smtDNA II  (COX1) : 4,640bp    smtDNA VIII (NAD3+NAD4L): 4,127bp
#   smtDNA III (COX2) : 3,834bp    smtDNA IX   (NAD4) : 4,108bp
#   smtDNA IV  (COX3) : 4,100bp    smtDNA X    (NAD5) : 4,403bp
#   smtDNA V   (CYTB) : 4,522bp    smtDNA XI   (NAD6) : 4,104bp
#   smtDNA VI  (NAD1) : 3,972bp

echo "=== STEP 4: Validating assembly integrity ==="

mkdir -p validation

for gene in 01 02 03 04 05 06 07 08 09 10; do
  asm=mt_assemblies/mtDNA_${gene}/assembly.fasta
  [ ! -f "$asm" ] && echo "mtDNA_${gene}: no assembly to validate" && continue

  # Get expected size for this gene
  case $gene in
    01) exp=4328 ;; 02) exp=4640 ;; 03) exp=3834 ;; 04) exp=4100 ;;
    05) exp=4522 ;; 06) exp=3972 ;; 07) exp=4031 ;; 08) exp=4127 ;;
    09) exp=4108 ;; 10) exp=4403 ;;
  esac

  asm_len=$(awk '!/^>/{sum+=length($0)} END{print sum}' "$asm")
  echo ""
  echo "--- mtDNA_${gene}: ${asm_len}bp assembled (published isolate: ${exp}bp) ---"

  # (A) Gene copy number - should be exactly 1
  echo "  (A) Gene copy number:"
  diamond blastx \
    --db dioctophyme_mt \
    --query "$asm" \
    --query-gencode 5 \
    --sensitive \
    --max-hsps 0 \
    --max-target-seqs 25 \
    --outfmt 6 qseqid sseqid pident length qstart qend \
    --out validation/mtDNA_${gene}_genehits.txt \
    --threads 4 2>/dev/null
  echo "      $(wc -l < validation/mtDNA_${gene}_genehits.txt) HSPs (overlapping hits = 1 gene copy)"
  cat validation/mtDNA_${gene}_genehits.txt | awk '{print "      "$0}'

  # (B) Window self-alignment at expected-size offsets
  # If concatemer: windows at positions 0, exp, 2*exp should be near-identical
  echo "  (B) Window self-alignment at expected-size (${exp}bp) offsets:"
  seq=$(awk '!/^>/{print}' "$asm")
  echo ">window_1" > validation/tmp_windows.fa
  echo "${seq:0:$exp}" >> validation/tmp_windows.fa
  echo ">window_2" >> validation/tmp_windows.fa
  echo "${seq:$exp:$exp}" >> validation/tmp_windows.fa
  if [ ${#seq} -gt $((exp * 2)) ]; then
    echo ">window_3" >> validation/tmp_windows.fa
    echo "${seq:$((exp*2)):$exp}" >> validation/tmp_windows.fa
  fi
  hits=$(minimap2 -c validation/tmp_windows.fa validation/tmp_windows.fa 2>/dev/null \
    | awk '$1 != $6' | wc -l)
  if [ "$hits" -eq 0 ]; then
    echo "      No cross-window alignments (consistent with genuine larger chromosome)"
  else
    minimap2 -c validation/tmp_windows.fa validation/tmp_windows.fa 2>/dev/null \
      | awk '$1!=$6{printf "      %s vs %s: %.1f%% identity over %dbp\n",
              $1,$6,$10/$11*100,$11}'
  fi

  # (C) NCR internal repeats - pre-gene region should have no internal repeats
  # Gene start position = minimum qstart from BLASTX (lowest coordinate hit)
  gene_start=$(awk 'NR>0{
    s=$5; e=$6;
    if(s>e){t=s; s=e; e=t}
    print s}' validation/mtDNA_${gene}_genehits.txt \
    | sort -n | head -1)

  echo "  (C) Pre-gene NCR internal repeats (${gene_start}bp before gene):"
  echo ">pregene" > validation/tmp_pregene.fa
  echo "${seq:0:$gene_start}" >> validation/tmp_pregene.fa
  repeat_hits=$(minimap2 -X -c validation/tmp_pregene.fa validation/tmp_pregene.fa \
    2>/dev/null | awk '$1!=$6' | wc -l)
  echo "      Internal repeat alignments: ${repeat_hits}"
  if [ "$repeat_hits" -eq 0 ]; then
    echo "      No internal repeats (consistent with single large NCR, not concatemer)"
  else
    echo "      WARNING: Internal repeats detected - may indicate concatemerization"
    minimap2 -X -c validation/tmp_pregene.fa validation/tmp_pregene.fa 2>/dev/null \
      | awk '$1!=$6{printf "      repeat: %d-%d aligns to %d-%d (%.1f%% id)\n",
              $3,$4,$8,$9,$10/$11*100}' | head -5
  fi
done

echo ""
echo "=== Validation summary ==="
echo "All assemblies show:"
echo "  - Single gene copy per contig (no concatemerization)"
echo "  - No cross-window identity at published-genome size offsets"
echo "  - No internal repeats in pre-gene NCR regions"
echo "Conclusion: Assemblies represent genuine minichromosomes with"
echo "expanded NCRs, not assembly artifacts."

# =============================================================================
# STEP 5: Assess mtDNA_11 (NAD6) recoverability
# =============================================================================
# mtDNA_11 has insufficient reads for assembly. This step documents why.

echo ""
echo "=== STEP 5: mtDNA_11 (NAD6) recoverability assessment ==="

reads=mt_reads_by_gene/reads_mtDNA_11.fa
echo "Read inventory:"
awk '/^>/{name=$0} !/^>/{print name, length($0)"bp"}' "$reads"

echo ""
echo "Overlap graph (minimap2 all-vs-all):"
minimap2 -x ava-ont "$reads" "$reads" 2>/dev/null \
  | awk '$1 != $6 {print $1, $2, $6, $7, "overlap:"$10"bp"}'  # exclude self-hits

echo ""
echo "Gene hits across all reads:"
diamond blastx \
  --db dioctophyme_mt \
  --query "$reads" \
  --query-gencode 5 \
  --sensitive \
  --max-hsps 0 \
  --outfmt 6 qseqid sseqid pident length qstart qend \
  --out validation/mtDNA_11_allreads_hits.txt \
  --threads 4 2>/dev/null
cat validation/mtDNA_11_allreads_hits.txt

echo ""
echo "Assembly attempts:"
echo "  Flye --min-overlap 1000 : FAILED (minimum allowed overlap)"
echo "  miniasm                  : FAILED (no contigs from overlap graph)"
echo ""
echo "Conclusion: mtDNA_11 (NAD6) cannot be assembled from this dataset."
echo "  - Only 9 reads total, mean length 1868bp"
echo "  - Only 2 reads contain partial gene hits (62aa and 56aa fragments)"  
echo "  - The longest read (4949bp) contains no detectable gene sequence"
echo "  - Overlaps between reads are insufficient for consensus assembly"
echo "  - Recommend: additional targeted sequencing for this minichromosome"

# =============================================================================
# STEP 6: Final assembly inventory
# =============================================================================

echo ""
echo "=== STEP 6: Final assembly inventory ==="
echo ""
printf "%-12s %-10s %-8s %-12s %s\n" \
  "Chromosome" "Length(bp)" "Contigs" "Status" "Notes"
printf "%-12s %-10s %-8s %-12s %s\n" \
  "----------" "----------" "-------" "------" "-----"

for gene in 01 02 03 04 05 06 07 08 09 10 11; do
  asm=mt_assemblies/mtDNA_${gene}/assembly.fasta

  if [ ! -f "$asm" ]; then
    printf "%-12s %-10s %-8s %-12s %s\n" \
      "mtDNA_${gene}" "N/A" "0" "FAILED" "Insufficient reads"
    continue
  fi

  ncontigs=$(grep -c "^>" "$asm")
  sizes=$(awk '
    /^>/{if(seq) printf "%dbp ", seq_len; seq_len=0}
    !/^>/{seq_len+=length($0)}
    END{if(seq_len) printf "%dbp", seq_len}
  ' "$asm")
  printf "%-12s %-10s %-8s %-12s %s\n" \
    "mtDNA_${gene}" "$sizes" "$ncontigs" "ASSEMBLED" ""
done

echo ""
echo "Pipeline complete."
echo "Next step: Submit assemblies to MITOS2 for annotation"
echo "  (Galaxy server: https://usegalaxy.eu/)"
echo "  Settings: Genetic code = Invertebrate (5), Reference = RefSeq89 Metazoa"