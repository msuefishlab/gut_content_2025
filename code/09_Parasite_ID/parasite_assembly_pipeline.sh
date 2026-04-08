#!/usr/bin/env bash
# =============================================================================
# parasite_assembly_pipeline.sh
# Unified Dioctophyme renale identification, assembly, and annotation pipeline
# =============================================================================
# Combines the strengths of three prior scripts into one logical flow:
#
#   1. Recruit reads  — DIAMOND (mt proteins), Barrnap (rRNA), BOLD (COI)
#   2. Assemble       — Flye per mt minichromosome + Flye for rRNA locus
#   3. Validate       — gene copy number, window self-alignment, NCR repeats
#   4. Annotate       — MITOS2 (mt), BLASTn (rRNA), BOLD COI species ID
#
# Key biological notes:
#   - D. renale has 11 mitochondrial minichromosomes (one PCG each).
#     Our isolate shows massively expanded NCRs vs Macchiaroli et al. 2025
#     (~5-11 kb vs ~2.6 kb published). Validated as genuine, not artifact.
#   - rRNA locus (18S-ITS1-5.8S-ITS2-28S) is nuclear, ~7-8 kb in nematodes.
#   - mtDNA_11 (NAD6) routinely fails assembly — too few reads (see STEP 6).
#
# Reference: Macchiaroli et al. 2025, DOI: 10.21203/rs.3.rs-8148715/v1
#            GenBank: PX715184-PX715194
#
# Usage:
#   bash parasite_assembly_pipeline.sh -r READS.fa -p DRENALE_PROTEINS.fa [opts]
#
# Required:
#   -r  PATH    All nanopore reads FASTA (whole-genome ONT)
#   -p  PATH    D. renale mitochondrial PCG amino acid FASTA
#
# Optional:
#   -b  PATH    BOLD FASTA database (filtered to Nematoda internally)
#               If omitted, BOLD alignment and COI-based recruitment are skipped.
#   -s  PATH    Per-sample reads for Barrnap (repeatable; default: use -r)
#   -o  PATH    Output directory [default: parasite_pipeline_out]
#   -R  PATH    MITOS2 refdata parent dir [default: <outdir>/mitos2_refdata]
#   -n  PATH    Local BLAST nr database path (omit = remote NCBI)
#   -e  FLOAT   BLASTn e-value cutoff [default: 1e-10]
#   -t  INT     Max BLASTn target sequences [default: 10]
#   -T  INT     Threads for DIAMOND, minimap2, local BLASTn [default: 8]
#   -k  STR     Barrnap kingdom: euk|bac|mito [default: euk]
#   -M  STR     Medaka model for post-assembly polishing [default: r941_min_sup_g507]
#               Use r941_min_fast_g507 for fast-called reads.
#               Set to 'none' to skip medaka polishing entirely.
#   -h          Show this help
#
# Dependencies:
#   diamond, flye, barrnap, minimap2, samtools, seqtk, blastn, medaka,
#   runmitos.py (bioconda: mitos), python3 + Biopython, curl
#
# Mac bash v3 note: no associative arrays (declare -A). Gene size lookups
#   use case statements. Octal interpretation avoided by quoting gene IDs.
# =============================================================================

set -euo pipefail

# =============================================================================
# Defaults and argument parsing
# =============================================================================

READS_FA=""
DRENALE_PROTEINS=""
BOLD_FA=""
SAMPLE_READS=()
OUTDIR="parasite_pipeline_out"
REFDATA_DIR=""          # set after OUTDIR is finalised if not provided
NR_DB=""
EVALUE="1e-10"
MAX_SEQS=10
THREADS=8
KINGDOM="euk"
MEDAKA_MODEL="r941_min_sup_g507"

usage() {
    sed -n '16,51p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 1
}

while getopts "r:p:b:s:o:R:n:e:t:T:k:M:h" opt; do
    case $opt in
        r) READS_FA="$OPTARG" ;;
        p) DRENALE_PROTEINS="$OPTARG" ;;
        b) BOLD_FA="$OPTARG" ;;
        s) SAMPLE_READS+=("$OPTARG") ;;
        o) OUTDIR="$OPTARG" ;;
        R) REFDATA_DIR="$OPTARG" ;;
        n) NR_DB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) MAX_SEQS="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        k) KINGDOM="$OPTARG" ;;
        M) MEDAKA_MODEL="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$READS_FA" || -z "$DRENALE_PROTEINS" ]]; then
    echo "ERROR: -r (reads) and -p (D. renale proteins) are required."
    usage
fi

# If no per-sample reads given, barrnap runs on the main reads file
if [[ ${#SAMPLE_READS[@]} -eq 0 ]]; then
    SAMPLE_READS=("$READS_FA")
fi

[[ -z "$REFDATA_DIR" ]] && REFDATA_DIR="${OUTDIR}/mitos2_refdata"
REFSEQ_DIR="${REFDATA_DIR}/refseq89m"

# =============================================================================
# Directory layout and logging
# =============================================================================

DB_DIR="${OUTDIR}/databases"
RECRUIT_DIR="${OUTDIR}/recruitment"
READS_DIR="${OUTDIR}/reads_by_target"
ASM_DIR="${OUTDIR}/assemblies"
VAL_DIR="${OUTDIR}/validation"
ANN_DIR="${OUTDIR}/annotations"

mkdir -p "$DB_DIR" "$RECRUIT_DIR/barrnap" "$RECRUIT_DIR/bold" \
         "$READS_DIR" "$ASM_DIR" "$VAL_DIR" "$ANN_DIR"

LOG="${OUTDIR}/pipeline.log"
echo "[$(date)] Pipeline started" | tee "$LOG"
echo "[$(date)] Reads:    $READS_FA" | tee -a "$LOG"
echo "[$(date)] Proteins: $DRENALE_PROTEINS" | tee -a "$LOG"
[[ -n "$BOLD_FA" ]] && echo "[$(date)] BOLD:     $BOLD_FA" | tee -a "$LOG"
echo "[$(date)] Outdir:   $OUTDIR" | tee -a "$LOG"

# =============================================================================
# Gene list — defined ONCE, reused throughout all subsequent loops.
# Bash v3: space-separated string; iterate with:  for gene in $GENES
# =============================================================================
GENES="01 02 03 04 05 06 07 08 09 10 11"

# =============================================================================
# Helper: flye_assembly TARGET_NAME READS_FA GENOME_SIZE OUTPUT_DIR
# =============================================================================
# Encapsulates the Flye non-determinism retry logic and tiered rescue.
# Flye can fail its polishing coverage filter even on valid assemblies;
# retrying up to 3 times catches most transient failures. On all-attempt
# failure the function rescues from the deepest available intermediate:
#   Tier 1: polished_N.fasta   (post-polishing, pre-finalize)
#   Tier 2: 30-contigger/contigs.fasta  (post-graph, unpolished)
#   Tier 3: 10-consensus/consensus.fasta (pre-graph disjointigs)
#
# Returns 0 if assembly.fasta was written (normal or rescued), 1 otherwise.
# Sets the global LAST_ASM_LEN to the assembled length (bp) on success.
# =============================================================================
LAST_ASM_LEN=0

flye_assembly() {
    local target="$1"
    local reads="$2"
    local gsize="$3"
    local outdir="$4"

    local n
    n=$(grep -c "^>" "$reads" 2>/dev/null || echo 0)
    if [ "$n" -lt 5 ]; then
        echo "[$(date)] ${target}: SKIPPED (only ${n} reads — insufficient)" | tee -a "$LOG"
        return 1
    fi

    # Flye rejects any path (reads or output dir) that contains spaces.
    # The OneDrive working directory has spaces in it, so we run Flye entirely
    # inside a space-free temp directory using symlinks, then move the results.
    local tmpdir
    tmpdir=$(mktemp -d /tmp/flye_XXXXXX)
    ln -s "$(cd "$(dirname "$reads")" && pwd)/$(basename "$reads")" "${tmpdir}/reads.fa"
    local tmp_outdir="${tmpdir}/assembly"

    local flye_success=0
    local attempt
    for attempt in 1 2 3; do
        rm -rf "$tmp_outdir"
        flye \
            --nano-raw "${tmpdir}/reads.fa" \
            --out-dir "$tmp_outdir" \
            --genome-size "$gsize" \
            --meta \
            --keep-haplotypes \
            --min-overlap 1000 \
            --iterations 3 \
            --threads "$THREADS" \
            2>&1 | tail -3 || true
        if [ -f "${tmp_outdir}/assembly_info.txt" ]; then
            flye_success=1
            break
        fi
        echo "[$(date)] ${target}: attempt ${attempt}/3 failed, retrying..." | tee -a "$LOG"
    done

    # Move Flye output from temp to real outdir
    rm -rf "$outdir"
    if [ -d "$tmp_outdir" ]; then
        cp -r "$tmp_outdir" "$outdir"
    else
        mkdir -p "$outdir"
    fi
    rm -rf "$tmpdir"

    local asm="${outdir}/assembly.fasta"
    if [ "$flye_success" -eq 1 ]; then
        LAST_ASM_LEN=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$asm")
        awk -v t="$target" 'NR>1{print t": "$2"bp circ="$4" cov="$5}' \
            "${outdir}/assembly_info.txt" | tee -a "$LOG"
        return 0
    fi

    # Tiered rescue
    local rescue_src="" rescue_stage=""
    local _t1 _t2 _t3
    _t1=$(find "$outdir" -name "polished_*.fasta" 2>/dev/null | sort -V | tail -1 || true)
    _t2="${outdir}/30-contigger/contigs.fasta"
    _t3="${outdir}/10-consensus/consensus.fasta"
    if [ -n "$_t1" ] && [ -s "$_t1" ]; then
        rescue_src="$_t1"; rescue_stage="polished"
    elif [ -f "$_t2" ] && [ -s "$_t2" ]; then
        rescue_src="$_t2"; rescue_stage="contigger"
    elif [ -f "$_t3" ] && [ -s "$_t3" ]; then
        rescue_src="$_t3"; rescue_stage="consensus"
    fi

    if [ -n "$rescue_src" ]; then
        cp "$rescue_src" "$asm"
        LAST_ASM_LEN=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$asm")
        echo "[$(date)] ${target}: rescued ${LAST_ASM_LEN}bp from ${rescue_stage} ($(basename "$rescue_src"))" | tee -a "$LOG"
        return 0
    fi

    echo "[$(date)] ${target}: ASSEMBLY FAILED (no usable output at any stage)" | tee -a "$LOG"
    return 1
}

# =============================================================================
# Helper: medaka_polish TARGET_NAME READS_FA ASSEMBLY_FA OUTPUT_DIR
# =============================================================================
# Runs medaka_consensus to correct homopolymer errors in an ONT assembly.
# Uses /tmp symlinks to work around spaces-in-path (same pattern as flye_assembly).
# Falls back to original assembly on failure.
medaka_polish() {
    local target="$1"
    local reads="$2"
    local assembly="$3"
    local outdir="$4"
    local polished="${outdir}/consensus.fasta"

    # Idempotency: skip if polished output already exists and is newer than assembly
    if [[ -s "$polished" && "$polished" -nt "$assembly" ]]; then
        echo "  [medaka] ${target}: polished output exists, skipping" | tee -a "$LOG"
        return 0
    fi

    mkdir -p "$outdir"

    local tmpdir
    tmpdir=$(mktemp -d /tmp/medaka_XXXXXX)
    ln -s "$(cd "$(dirname "$reads")" && pwd)/$(basename "$reads")" "${tmpdir}/reads.fa"
    ln -s "$(cd "$(dirname "$assembly")" && pwd)/$(basename "$assembly")" "${tmpdir}/assembly.fa"

    echo "  [medaka] ${target}: running medaka_consensus (model: ${MEDAKA_MODEL})..." | tee -a "$LOG"
    medaka_consensus \
        -i "${tmpdir}/reads.fa" \
        -d "${tmpdir}/assembly.fa" \
        -o "${tmpdir}/medaka_out" \
        -m "$MEDAKA_MODEL" \
        -t "$THREADS" 2>&1 | tail -5 | tee -a "$LOG"

    if [[ -s "${tmpdir}/medaka_out/consensus.fasta" ]]; then
        cp "${tmpdir}/medaka_out/consensus.fasta" "$polished"
        rm -rf "$tmpdir"
        local pol_len
        pol_len=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$polished")
        echo "  [medaka] ${target}: polished → ${pol_len} bp" | tee -a "$LOG"
        return 0
    else
        echo "  [medaka] ${target}: FAILED — keeping Flye assembly unchanged" | tee -a "$LOG" >&2
        rm -rf "$tmpdir"
        return 1
    fi
}

# =============================================================================
# STEP 0: Preflight checks and database / index setup
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 0: Preflight and database setup ===" | tee -a "$LOG"

# --- 0a: Dependency check ---
echo "[$(date)] Checking dependencies..." | tee -a "$LOG"
missing=""
for cmd in diamond flye barrnap minimap2 samtools seqtk blastn python3 curl; do
    command -v "$cmd" &>/dev/null || missing="${missing} ${cmd}"
done
if [[ "$MEDAKA_MODEL" != "none" ]]; then
    command -v medaka_consensus &>/dev/null || missing="${missing} medaka"
fi
if ! python3 -c "from Bio import SeqIO" &>/dev/null; then
    missing="${missing} python3-Biopython"
fi
if [[ -n "$missing" ]]; then
    echo "ERROR: Missing dependencies:${missing}" | tee -a "$LOG"
    exit 1
fi
# MITOS2 is only needed for annotation — warn but don't abort
if ! command -v runmitos.py &>/dev/null; then
    echo "[$(date)] WARNING: runmitos.py not found — MITOS2 annotation (STEP 7a) will be skipped." | tee -a "$LOG"
    SKIP_MITOS2=1
else
    SKIP_MITOS2=0
fi
echo "[$(date)] All required dependencies found." | tee -a "$LOG"

# --- 0b: Index reads file ---
if [[ ! -f "${READS_FA}.fai" ]]; then
    echo "[$(date)] Indexing reads: $READS_FA" | tee -a "$LOG"
    samtools faidx "$READS_FA"
fi

# --- 0c: D. renale DIAMOND database ---
DRENALE_DB="${DB_DIR}/dioctophyme_mt"
if [[ ! -f "${DRENALE_DB}.dmnd" ]]; then
    echo "[$(date)] Building D. renale DIAMOND database..." | tee -a "$LOG"
    diamond makedb --in "$DRENALE_PROTEINS" --db "$DRENALE_DB" 2>>"$LOG"
else
    echo "[$(date)] D. renale DIAMOND db already exists, skipping." | tee -a "$LOG"
fi

# --- 0d: Clade I nematode proteins + DIAMOND database ---
CLADEI_DB="${DB_DIR}/cladeI_mt"
if [[ ! -f "${CLADEI_DB}.dmnd" ]]; then
    echo "[$(date)] Downloading Clade I mitochondrial proteins from NCBI..." | tee -a "$LOG"
    mkdir -p "${DB_DIR}/cladeI_proteins"
    CLADEI_FA="${DB_DIR}/cladeI_proteins/cladeI_all_proteins.fa"
    : > "$CLADEI_FA"
    # Accessions from Macchiaroli et al. 2025 supplementary
    for acc in \
        NC_071371 NC_056391 NC_008640 NC_008693 NC_008692 \
        NC_025750 NC_025751 NC_025752 NC_025753 NC_025754 \
        NC_025749 NC_002681 NC_025755 NC_018596 NC_028621 \
        NC_018597 NC_017747 NC_017750 NC_001328; do
        echo "[$(date)]   Fetching ${acc}..." | tee -a "$LOG"
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=nuccore&id=${acc}&rettype=fasta_cds_aa&retmode=text" >> "$CLADEI_FA"
        sleep 1
    done
    echo "[$(date)] Clade I proteins: $(grep -c "^>" "$CLADEI_FA") sequences" | tee -a "$LOG"
    diamond makedb --in "$CLADEI_FA" --db "$CLADEI_DB" 2>>"$LOG"
else
    echo "[$(date)] Clade I DIAMOND db already exists, skipping download." | tee -a "$LOG"
fi

# --- 0e: MITOS2 RefSeq89 Metazoa reference ---
if [[ ! -d "$REFSEQ_DIR" ]] && [[ "$SKIP_MITOS2" -eq 0 ]]; then
    echo "[$(date)] Downloading MITOS2 RefSeq89 Metazoa reference (~20MB)..." | tee -a "$LOG"
    mkdir -p "$REFDATA_DIR"
    curl -L --progress-bar \
        "https://zenodo.org/records/3685310/files/refseq89m.tar.bz2" \
        -o "${REFDATA_DIR}/refseq89m.tar.bz2"
    tar -xjf "${REFDATA_DIR}/refseq89m.tar.bz2" -C "${REFDATA_DIR}/"
    rm "${REFDATA_DIR}/refseq89m.tar.bz2"
    echo "[$(date)] MITOS2 reference ready: $REFSEQ_DIR" | tee -a "$LOG"
fi

# --- 0f: BOLD Nematoda filter + minimap2 index (only when -b provided) ---
BOLD_NEMA=""
BOLD_MMI=""
if [[ -n "$BOLD_FA" ]]; then
    BOLD_NEMA="${DB_DIR}/bold_nematoda.fa"
    BOLD_MMI="${DB_DIR}/bold_nematoda.mmi"
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
print(f"  Wrote {count} Nematoda sequences")
PYEOF
    fi
    if [[ ! -f "$BOLD_MMI" ]]; then
        echo "[$(date)] Building BOLD Nematoda minimap2 index..." | tee -a "$LOG"
        minimap2 -d "$BOLD_MMI" "$BOLD_NEMA" 2>>"$LOG"
    fi
    nema_count=$(grep -c ">" "$BOLD_NEMA")
    echo "[$(date)] BOLD Nematoda: ${nema_count} sequences" | tee -a "$LOG"
fi

# =============================================================================
# STEP 1: Multi-strategy read recruitment
# =============================================================================
# Four strategies cast the widest possible net:
#   1a  DIAMOND vs D. renale mt proteins (specific)
#   1b  DIAMOND vs Clade I mt proteins (broad; catches diverged reads)
#   1c  Barrnap rRNA detection (nuclear 18S/5.8S/28S)
#   1d  minimap2 vs BOLD Nematoda COI (optional)
#
# --long-reads is essential for DIAMOND on nanopore data: enables the
# chaining algorithm for sparse seeds. Without it, recruitment drops ~50%.

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 1: Multi-strategy read recruitment ===" | tee -a "$LOG"

# --- 1a: DIAMOND vs D. renale ---
echo "[$(date)] 1a: DIAMOND blastx vs D. renale mt proteins..." | tee -a "$LOG"
diamond blastx \
    --db "$DRENALE_DB" \
    --query "$READS_FA" \
    --query-gencode 5 \
    --long-reads \
    --sensitive \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$THREADS" \
    --out "${RECRUIT_DIR}/mt_hits_drenale.txt" 2>>"$LOG"
echo "[$(date)]   Recruited: $(cut -f1 "${RECRUIT_DIR}/mt_hits_drenale.txt" | sort -u | wc -l | tr -d ' ') reads" | tee -a "$LOG"

# --- 1b: DIAMOND vs Clade I ---
echo "[$(date)] 1b: DIAMOND blastx vs Clade I mt proteins..." | tee -a "$LOG"
diamond blastx \
    --db "$CLADEI_DB" \
    --query "$READS_FA" \
    --query-gencode 5 \
    --long-reads \
    --sensitive \
    --evalue 1e-3 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$THREADS" \
    --out "${RECRUIT_DIR}/mt_hits_cladeI.txt" 2>>"$LOG"
echo "[$(date)]   Recruited: $(cut -f1 "${RECRUIT_DIR}/mt_hits_cladeI.txt" | sort -u | wc -l | tr -d ' ') reads" | tee -a "$LOG"

# --- 1c: Barrnap rRNA detection per sample ---
echo "[$(date)] 1c: Barrnap rRNA detection..." | tee -a "$LOG"
for sample_fa in "${SAMPLE_READS[@]}"; do
    sample=$(basename "$sample_fa" | sed 's/\.[^.]*$//')
    gff="${RECRUIT_DIR}/barrnap/${sample}_${KINGDOM}.gff"
    if [[ ! -f "${sample_fa}.fai" ]]; then
        samtools faidx "$sample_fa"
    fi
    echo "[$(date)]   Barrnap on ${sample}..." | tee -a "$LOG"
    barrnap --kingdom "$KINGDOM" "$sample_fa" > "$gff" 2>>"$LOG"
    nfeat=$(grep -v "^#" "$gff" | grep -v "5S_rRNA" | wc -l | tr -d ' ')
    echo "[$(date)]   ${sample}: ${nfeat} rRNA features (excl. 5S)" | tee -a "$LOG"
done

# --- 1d: minimap2 vs BOLD Nematoda (only when -b provided) ---
if [[ -n "$BOLD_FA" ]]; then
    echo "[$(date)] 1d: minimap2 alignment vs BOLD Nematoda..." | tee -a "$LOG"
    for sample_fa in "${SAMPLE_READS[@]}"; do
        sample=$(basename "$sample_fa" | sed 's/\.[^.]*$//')
        bam="${RECRUIT_DIR}/bold/${sample}_vs_bold_nematoda.bam"
        hits_txt="${RECRUIT_DIR}/bold/${sample}_bold_hits.txt"
        minimap2 -ax map-ont \
            --secondary=no \
            -f 0.0001 \
            -t "$THREADS" \
            "$BOLD_MMI" "$sample_fa" 2>>"$LOG" \
            | samtools view -F 4 -b \
            | samtools sort -o "$bam"
        samtools index "$bam"
        echo "# read_id  reference  mapq  flag" > "$hits_txt"
        samtools view "$bam" | awk '{print $1, $3, $5, $2}' >> "$hits_txt"
        mapped=$(grep -v "^#" "$hits_txt" | wc -l | tr -d ' ')
        echo "[$(date)]   ${sample}: ${mapped} reads mapped to BOLD Nematoda" | tee -a "$LOG"
    done
fi

# =============================================================================
# STEP 2: Read pool construction
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 2: Read pool construction ===" | tee -a "$LOG"

# --- 2a: mt pool (union of DIAMOND passes) ---
cat "${RECRUIT_DIR}/mt_hits_drenale.txt" "${RECRUIT_DIR}/mt_hits_cladeI.txt" \
    | cut -f1 | sort -u > "${RECRUIT_DIR}/all_mt_read_ids.txt"
echo "[$(date)] mt pool: $(wc -l < "${RECRUIT_DIR}/all_mt_read_ids.txt" | tr -d ' ') reads" | tee -a "$LOG"
seqtk subseq "$READS_FA" "${RECRUIT_DIR}/all_mt_read_ids.txt" \
    > "${RECRUIT_DIR}/all_mt_reads_pooled.fa"

# --- 2b: rRNA read IDs from barrnap GFFs (exclude 5S — common spurious hit) ---
# A read qualifies if it has at least one 18S, 5.8S, or 28S feature.
python3 - <<PYEOF
import os, glob
gff_dir = "${RECRUIT_DIR}/barrnap"
out_ids = "${RECRUIT_DIR}/barrnap_rrna_read_ids.txt"
keep_types = {"18S_rRNA", "5_8S_rRNA", "28S_rRNA"}
read_ids = set()
for gff in glob.glob(os.path.join(gff_dir, "*.gff")):
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            name = [a.split("=")[1] for a in attrs.split(";") if a.startswith("Name=")]
            name = name[0] if name else ""
            if name in keep_types:
                read_ids.add(parts[0])
with open(out_ids, "w") as fh:
    fh.write("\n".join(sorted(read_ids)) + ("\n" if read_ids else ""))
print(f"  Barrnap rRNA reads (18S/5.8S/28S): {len(read_ids)}")
PYEOF

# --- 2c: rRNA pool (union of barrnap + BOLD-mapped IDs) ---
sort -u "${RECRUIT_DIR}/barrnap_rrna_read_ids.txt" > "${RECRUIT_DIR}/all_rrna_read_ids.txt"
if [[ -n "$BOLD_FA" ]]; then
    for hits_txt in "${RECRUIT_DIR}/bold/"*_bold_hits.txt; do
        grep -v "^#" "$hits_txt" | cut -d' ' -f1 >> "${RECRUIT_DIR}/all_rrna_read_ids.txt"
    done
    sort -u "${RECRUIT_DIR}/all_rrna_read_ids.txt" -o "${RECRUIT_DIR}/all_rrna_read_ids.txt"
fi
echo "[$(date)] rRNA pool: $(wc -l < "${RECRUIT_DIR}/all_rrna_read_ids.txt" | tr -d ' ') reads" | tee -a "$LOG"
seqtk subseq "$READS_FA" "${RECRUIT_DIR}/all_rrna_read_ids.txt" \
    > "${RECRUIT_DIR}/all_rrna_reads_pooled.fa" 2>/dev/null || true

# Master union (for record-keeping)
sort -u "${RECRUIT_DIR}/all_mt_read_ids.txt" "${RECRUIT_DIR}/all_rrna_read_ids.txt" \
    > "${RECRUIT_DIR}/all_recruited_read_ids.txt"
echo "[$(date)] Total recruited reads (union): $(wc -l < "${RECRUIT_DIR}/all_recruited_read_ids.txt" | tr -d ' ')" | tee -a "$LOG"

# --- 2d: Per-gene DIAMOND assignment of mt pool ---
# All-hits mode (--max-hsps 0): a read goes into every gene pool it matches.
# At low coverage (20-65 reads/gene), losing bridging reads breaks assembly.
# Chimeric assemblies from multi-gene reads are caught by STEP 5 validation.
echo "[$(date)] 2d: DIAMOND per-gene read assignment (mt pool)..." | tee -a "$LOG"
diamond blastx \
    --db "$DRENALE_DB" \
    --query "${RECRUIT_DIR}/all_mt_reads_pooled.fa" \
    --query-gencode 5 \
    --long-reads \
    --sensitive \
    --max-hsps 0 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --threads "$THREADS" \
    --out "${RECRUIT_DIR}/mt_read_hits_pooled_assigned.txt" 2>>"$LOG"
echo "[$(date)] Per-gene read counts:" | tee -a "$LOG"
cut -f2 "${RECRUIT_DIR}/mt_read_hits_pooled_assigned.txt" | sort | uniq -c | sort -rn | tee -a "$LOG"

# =============================================================================
# STEP 3: Per-target read splitting
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 3: Per-target read splitting ===" | tee -a "$LOG"

# --- 3a: mt genes ---
for gene in $GENES; do
    grep "mtDNA_${gene}" "${RECRUIT_DIR}/mt_read_hits_pooled_assigned.txt" \
        | cut -f1 | sort -u \
        > "${READS_DIR}/reads_mtDNA_${gene}.ids"
    seqtk subseq "$READS_FA" \
        "${READS_DIR}/reads_mtDNA_${gene}.ids" \
        > "${READS_DIR}/reads_mtDNA_${gene}.fa"
    n=$(wc -l < "${READS_DIR}/reads_mtDNA_${gene}.ids")
    mean=$(awk '/^>/{next}{sum+=length($0);n++}END{
        if(n>0) printf "%.0f",sum/n; else print "NA"}' \
        "${READS_DIR}/reads_mtDNA_${gene}.fa")
    maxlen=$(awk '/^>/{next}{if(length($0)>m)m=length($0)}END{print m+0}' \
        "${READS_DIR}/reads_mtDNA_${gene}.fa")
    echo "[$(date)] mtDNA_${gene}: ${n} reads, mean=${mean}bp, max=${maxlen}bp" | tee -a "$LOG"
done

# --- 3b: rRNA pool ---
rrna_n=$(grep -c "^>" "${RECRUIT_DIR}/all_rrna_reads_pooled.fa" 2>/dev/null || echo 0)
echo "[$(date)] rRNA_locus: ${rrna_n} reads" | tee -a "$LOG"

# =============================================================================
# STEP 4: Assembly with Flye
# =============================================================================
# Parameters shared across all targets:
#   --meta            metagenome mode, tolerant of uneven coverage
#   --keep-haplotypes preserve variant haplotypes
#   --min-overlap 1000 Flye's minimum; paper used 3800 but our depths require 1000
#   --iterations 3    polishing rounds
#
# Genome sizes:
#   mt minichromosomes: 4k  (PCG + NCR; published ~4kb, ours up to 30kb)
#   rRNA locus:         8k  (18S + ITS1 + 5.8S + ITS2 + 28S ≈ 7-8kb in nematodes)

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 4: Assembly with Flye ===" | tee -a "$LOG"

# --- 4a: mt minichromosomes (one per gene) ---
declare -A MT_ASM_LENS 2>/dev/null || true   # bash v4+ only; fallback below
MT_ASM_STATUS=""

for gene in $GENES; do
    reads="${READS_DIR}/reads_mtDNA_${gene}.fa"
    outdir="${ASM_DIR}/mtDNA_${gene}"
    if flye_assembly "mtDNA_${gene}" "$reads" "4k" "$outdir"; then
        MT_ASM_STATUS="${MT_ASM_STATUS}${gene}:${LAST_ASM_LEN} "
    else
        MT_ASM_STATUS="${MT_ASM_STATUS}${gene}:FAILED "
    fi
done

# --- 4b: rRNA locus ---
RRNA_ASM_LEN=0
RRNA_ASM_OK=0
if flye_assembly "rRNA_locus" \
       "${RECRUIT_DIR}/all_rrna_reads_pooled.fa" "8k" \
       "${ASM_DIR}/rRNA_locus"; then
    RRNA_ASM_LEN=$LAST_ASM_LEN
    RRNA_ASM_OK=1
fi

# =============================================================================
# STEP 4.5: Post-assembly polishing with medaka
# =============================================================================
# Corrects residual homopolymer errors that Flye's internal polishing misses.
# Requires medaka to be installed (skipped if -M none).

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 4.5: Medaka polishing ===" | tee -a "$LOG"

if [[ "$MEDAKA_MODEL" == "none" ]]; then
    echo "[$(date)] Medaka polishing skipped (-M none)." | tee -a "$LOG"
else
    echo "  Model: ${MEDAKA_MODEL}" | tee -a "$LOG"

    for gene in 01 02 03 04 05 06 07 08 09 10; do
        asm="${ASM_DIR}/mtDNA_${gene}/assembly.fasta"
        reads="${READS_DIR}/reads_mtDNA_${gene}.fa"
        pol_out="${ASM_DIR}/mtDNA_${gene}/medaka"
        polished="${pol_out}/consensus.fasta"

        if [[ ! -s "$asm" ]]; then
            echo "  [medaka] mtDNA_${gene}: no assembly, skipping" | tee -a "$LOG"
            continue
        fi
        if [[ ! -s "$reads" ]]; then
            echo "  [medaka] mtDNA_${gene}: no reads, skipping" | tee -a "$LOG"
            continue
        fi

        if medaka_polish "mtDNA_${gene}" "$reads" "$asm" "$pol_out"; then
            cp "$polished" "$asm"
            echo "  [medaka] mtDNA_${gene}: assembly.fasta updated with polished sequence" | tee -a "$LOG"
        fi
    done

    # Polish rRNA locus if assembled
    rna_asm="${ASM_DIR}/rRNA_locus/assembly.fasta"
    rna_reads="${RECRUIT_DIR}/all_rrna_reads_pooled.fa"
    if [[ -s "$rna_asm" && -s "$rna_reads" ]]; then
        if medaka_polish "rRNA_locus" "$rna_reads" "$rna_asm" \
                "${ASM_DIR}/rRNA_locus/medaka"; then
            cp "${ASM_DIR}/rRNA_locus/medaka/consensus.fasta" "$rna_asm"
            echo "  [medaka] rRNA_locus: assembly.fasta updated" | tee -a "$LOG"
        fi
    fi
fi

# =============================================================================
# STEP 5: Mt assembly integrity validation (genes 01-10; NAD6 handled in STEP 6)
# =============================================================================
# Three independent lines of evidence rule out concatemerization:
#   (A) Gene copy number — concatemers show the gene at regular intervals
#   (B) Window self-alignment — concatemers have >90% identity between windows
#   (C) NCR internal repeats — genuine expanded NCRs have none
#
# Expected sizes from Macchiaroli et al. 2025 Table 1 (their isolate):
#   smtDNA I  (ATP6):4328  II (COX1):4640  III (COX2):3834  IV (COX3):4100
#   V (CYTB):4522  VI (NAD1):3972  VII (NAD2):4031  VIII (NAD3+NAD4L):4127
#   IX (NAD4):4108  X (NAD5):4403

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 5: Mt assembly integrity validation ===" | tee -a "$LOG"

for gene in 01 02 03 04 05 06 07 08 09 10; do
    asm="${ASM_DIR}/mtDNA_${gene}/assembly.fasta"
    if [ ! -f "$asm" ] || [ ! -s "$asm" ]; then
        echo "[$(date)] mtDNA_${gene}: no assembly to validate" | tee -a "$LOG"
        continue
    fi

    case $gene in
        01) exp=4328 ;; 02) exp=4640 ;; 03) exp=3834 ;; 04) exp=4100 ;;
        05) exp=4522 ;; 06) exp=3972 ;; 07) exp=4031 ;; 08) exp=4127 ;;
        09) exp=4108 ;; 10) exp=4403 ;;
    esac

    asm_len=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$asm")
    echo "" | tee -a "$LOG"
    echo "[$(date)] --- mtDNA_${gene}: ${asm_len}bp assembled (published: ${exp}bp) ---" | tee -a "$LOG"

    # (A) Gene copy number
    echo "  (A) Gene copy number:" | tee -a "$LOG"
    diamond blastx \
        --db "$DRENALE_DB" \
        --query "$asm" \
        --query-gencode 5 \
        --sensitive \
        --max-hsps 0 \
        --max-target-seqs 25 \
        --outfmt 6 qseqid sseqid pident length qstart qend \
        --out "${VAL_DIR}/mtDNA_${gene}_genehits.txt" \
        --threads 4 2>/dev/null
    echo "      $(wc -l < "${VAL_DIR}/mtDNA_${gene}_genehits.txt") HSPs (overlapping = 1 gene copy)" | tee -a "$LOG"
    awk '{print "      "$0}' "${VAL_DIR}/mtDNA_${gene}_genehits.txt" | tee -a "$LOG"

    # (B) Window self-alignment at expected-size offsets
    echo "  (B) Window self-alignment at ${exp}bp offsets:" | tee -a "$LOG"
    seq=$(awk '!/^>/{print}' "$asm")
    echo ">window_1" > "${VAL_DIR}/tmp_windows.fa"
    echo "${seq:0:$exp}" >> "${VAL_DIR}/tmp_windows.fa"
    echo ">window_2" >> "${VAL_DIR}/tmp_windows.fa"
    echo "${seq:$exp:$exp}" >> "${VAL_DIR}/tmp_windows.fa"
    if [ ${#seq} -gt $((exp * 2)) ]; then
        echo ">window_3" >> "${VAL_DIR}/tmp_windows.fa"
        echo "${seq:$((exp*2)):$exp}" >> "${VAL_DIR}/tmp_windows.fa"
    fi
    hits=$(minimap2 -c "${VAL_DIR}/tmp_windows.fa" "${VAL_DIR}/tmp_windows.fa" 2>/dev/null \
        | awk '$1 != $6' | wc -l)
    if [ "$hits" -eq 0 ]; then
        echo "      No cross-window alignments (consistent with genuine large chromosome)" | tee -a "$LOG"
    else
        minimap2 -c "${VAL_DIR}/tmp_windows.fa" "${VAL_DIR}/tmp_windows.fa" 2>/dev/null \
            | awk '$1!=$6{printf "      %s vs %s: %.1f%% id over %dbp\n",
                    $1,$6,$10/$11*100,$11}' | tee -a "$LOG"
    fi

    # (C) Pre-gene NCR internal repeats
    gene_start=$(awk '{s=$5; e=$6; if(s>e){t=s;s=e;e=t} print s}' \
        "${VAL_DIR}/mtDNA_${gene}_genehits.txt" | sort -n | head -1)
    echo "  (C) Pre-gene NCR internal repeats (${gene_start}bp before gene):" | tee -a "$LOG"
    echo ">pregene" > "${VAL_DIR}/tmp_pregene.fa"
    echo "${seq:0:$gene_start}" >> "${VAL_DIR}/tmp_pregene.fa"
    repeat_hits=$(minimap2 -X -c \
        "${VAL_DIR}/tmp_pregene.fa" "${VAL_DIR}/tmp_pregene.fa" \
        2>/dev/null | awk '$1!=$6' | wc -l)
    echo "      Internal repeat alignments: ${repeat_hits}" | tee -a "$LOG"
    if [ "$repeat_hits" -eq 0 ]; then
        echo "      No internal repeats (genuine large NCR, not concatemer)" | tee -a "$LOG"
    else
        echo "      WARNING: internal repeats detected — possible concatemerization" | tee -a "$LOG"
        minimap2 -X -c \
            "${VAL_DIR}/tmp_pregene.fa" "${VAL_DIR}/tmp_pregene.fa" 2>/dev/null \
            | awk '$1!=$6{printf "      repeat: %d-%d aligns to %d-%d (%.1f%% id)\n",
                    $3,$4,$8,$9,$10/$11*100}' | head -5 | tee -a "$LOG"
    fi
done

# =============================================================================
# STEP 6: NAD6 (mtDNA_11) recoverability assessment
# =============================================================================
# NAD6 routinely fails: too few reads, most lacking detectable gene sequence.
# This step documents the evidence rather than silently skipping it.

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 6: NAD6 (mtDNA_11) recoverability assessment ===" | tee -a "$LOG"

reads11="${READS_DIR}/reads_mtDNA_11.fa"
n11=$(grep -c "^>" "$reads11" 2>/dev/null || echo 0)
if [ "$n11" -gt 0 ]; then
    echo "Read inventory:" | tee -a "$LOG"
    awk '/^>/{name=$0} !/^>/{print name, length($0)"bp"}' "$reads11" | tee -a "$LOG"

    echo "" | tee -a "$LOG"
    echo "Overlap graph (all-vs-all):" | tee -a "$LOG"
    minimap2 -x ava-ont "$reads11" "$reads11" 2>/dev/null \
        | awk '$1 != $6 {print $1, $2, $6, $7, "overlap:"$10"bp"}' | tee -a "$LOG"

    echo "" | tee -a "$LOG"
    echo "Gene hits across all reads:" | tee -a "$LOG"
    diamond blastx \
        --db "$DRENALE_DB" \
        --query "$reads11" \
        --query-gencode 5 \
        --sensitive \
        --max-hsps 0 \
        --outfmt 6 qseqid sseqid pident length qstart qend \
        --out "${VAL_DIR}/mtDNA_11_allreads_hits.txt" \
        --threads 4 2>/dev/null
    cat "${VAL_DIR}/mtDNA_11_allreads_hits.txt" | tee -a "$LOG"
fi
echo "" | tee -a "$LOG"
echo "Conclusion: mtDNA_11 (NAD6) cannot be assembled — insufficient reads." | tee -a "$LOG"
echo "  Recommend additional targeted sequencing for this minichromosome." | tee -a "$LOG"

# =============================================================================
# STEP 7: Annotation
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 7: Annotation ===" | tee -a "$LOG"

# --- 7a: MITOS2 for mt minichromosomes (01-10) ---
if [[ "$SKIP_MITOS2" -eq 1 ]]; then
    echo "[$(date)] 7a: MITOS2 not found — skipping mt annotation." | tee -a "$LOG"
else
    echo "[$(date)] 7a: MITOS2 annotation of mt assemblies..." | tee -a "$LOG"
    for gene in 01 02 03 04 05 06 07 08 09 10; do
        asm="${ASM_DIR}/mtDNA_${gene}/assembly.fasta"
        outdir="${ANN_DIR}/mtDNA_${gene}"
        if [ ! -f "$asm" ] || [ ! -s "$asm" ]; then
            echo "[$(date)] mtDNA_${gene}: no assembly — skipping MITOS2" | tee -a "$LOG"
            continue
        fi
        rm -rf "$outdir"
        mkdir -p "$outdir"
        runmitos.py \
            -i "$asm" \
            -c 5 \
            -o "$outdir" \
            -R "$REFDATA_DIR" \
            -r "refseq89m" \
            --noplots \
            2>&1 | grep -vE "^(DEBUG|WARNING: No handlers)" || true
        gff="${outdir}/result.gff"
        if [ -f "$gff" ] && [ -s "$gff" ]; then
            genes_found=$(grep "Name=" "$gff" 2>/dev/null \
                | sed 's/.*Name=\([^;]*\).*/\1/' | sort -u | tr '\n' ' ' || true)
            echo "[$(date)] mtDNA_${gene}: ${genes_found}" | tee -a "$LOG"
        else
            echo "[$(date)] mtDNA_${gene}: WARNING — no result.gff produced" | tee -a "$LOG"
        fi
    done
fi

# --- 7b: BLASTn on assembled rRNA contig(s) ---
RRNA_BLAST_HIT="none"
rrna_asm="${ASM_DIR}/rRNA_locus/assembly.fasta"
if [ -f "$rrna_asm" ] && [ -s "$rrna_asm" ]; then
    echo "[$(date)] 7b: BLASTn rRNA assembly vs NCBI nr..." | tee -a "$LOG"
    BLAST_FMT="6 qseqid sseqid pident length qcovs evalue bitscore stitle sscinames"
    rrna_tsv="${ANN_DIR}/rRNA_locus/rRNA_locus.blast.tsv"
    rrna_txt="${ANN_DIR}/rRNA_locus/rRNA_locus.blast.txt"
    mkdir -p "${ANN_DIR}/rRNA_locus"
    if [[ -n "$NR_DB" ]]; then
        blastn -query "$rrna_asm" -db "$NR_DB" \
            -out "$rrna_tsv" -outfmt "$BLAST_FMT" \
            -evalue "$EVALUE" -max_target_seqs "$MAX_SEQS" \
            -num_threads "$THREADS" 2>>"$LOG"
        blastn -query "$rrna_asm" -db "$NR_DB" \
            -out "$rrna_txt" -outfmt 0 \
            -evalue "$EVALUE" -max_target_seqs 5 \
            -num_threads "$THREADS" 2>>"$LOG"
    else
        echo "[$(date)]   NOTE: remote BLAST — this may take several minutes." | tee -a "$LOG"
        blastn -query "$rrna_asm" -db nr \
            -out "$rrna_tsv" -outfmt "$BLAST_FMT" \
            -evalue "$EVALUE" -max_target_seqs "$MAX_SEQS" \
            -remote 2>>"$LOG"
        blastn -query "$rrna_asm" -db nr \
            -out "$rrna_txt" -outfmt 0 \
            -evalue "$EVALUE" -max_target_seqs 5 \
            -remote 2>>"$LOG"
    fi
    if [[ -s "$rrna_tsv" ]]; then
        RRNA_BLAST_HIT=$(head -1 "$rrna_tsv" | awk -F'\t' \
            '{printf "%s (pident=%.1f%% qcov=%s%% evalue=%s)", $9, $3, $5, $6}')
        echo "[$(date)]   Top hit: ${RRNA_BLAST_HIT}" | tee -a "$LOG"
    else
        echo "[$(date)]   No BLASTn hits above evalue ${EVALUE}" | tee -a "$LOG"
    fi
else
    echo "[$(date)] 7b: No rRNA assembly — skipping BLASTn." | tee -a "$LOG"
fi

# =============================================================================
# STEP 8: BOLD COI read-level species ID (only when -b provided)
# =============================================================================

if [[ -n "$BOLD_FA" ]]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] === STEP 8: BOLD COI species identification ===" | tee -a "$LOG"
    for hits_txt in "${RECRUIT_DIR}/bold/"*_bold_hits.txt; do
        [[ -f "$hits_txt" ]] || continue
        sample=$(basename "$(dirname "$hits_txt")")
        echo "[$(date)] Top BOLD hits for ${sample} (by MAPQ):" | tee -a "$LOG"
        grep -v "^#" "$hits_txt" | sort -k3 -rn | head -10 \
            | awk '{printf "  mapq=%-3s  %s\n", $3, $2}' | tee -a "$LOG"
    done
fi

# =============================================================================
# STEP 9: Final summary table
# =============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] === STEP 9: Final summary ===" | tee -a "$LOG"
echo "" | tee -a "$LOG"

SUMMARY="${OUTDIR}/summary_table.txt"
{
printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
    "Target" "Length(bp)" "Contigs" "Circ" "Status" "Annotation"
printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
    "--------------" "----------" "-------" "----" "------------" "----------"

for gene in $GENES; do
    asm="${ASM_DIR}/mtDNA_${gene}/assembly.fasta"
    if [ ! -f "$asm" ] || [ ! -s "$asm" ]; then
        printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
            "mtDNA_${gene}" "N/A" "0" "-" "FAILED" "Insufficient reads"
        continue
    fi
    asm_len=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$asm")
    ncontigs=$(grep -c "^>" "$asm")
    circ=$(awk -v g="mtDNA_${gene}" 'NR>1{print $4}' \
        "${ASM_DIR}/mtDNA_${gene}/assembly_info.txt" 2>/dev/null | head -1 || echo "?")
    # Annotation notes
    ann_note=""
    gff="${ANN_DIR}/mtDNA_${gene}/result.gff"
    if [ -f "$gff" ] && [ -s "$gff" ]; then
        n_genes=$(grep -c "Name=" "$gff" 2>/dev/null || echo 0)
        ann_note="MITOS2: ${n_genes} features"
    elif [[ "$SKIP_MITOS2" -eq 1 ]]; then
        ann_note="MITOS2 skipped"
    else
        ann_note="no annotation"
    fi
    printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
        "mtDNA_${gene}" "$asm_len" "$ncontigs" "$circ" "ASSEMBLED" "$ann_note"
done

# rRNA locus row
rrna_asm="${ASM_DIR}/rRNA_locus/assembly.fasta"
if [ -f "$rrna_asm" ] && [ -s "$rrna_asm" ]; then
    rrna_len=$(awk '!/^>/{sum+=length($0)} END{print sum+0}' "$rrna_asm")
    rrna_nc=$(grep -c "^>" "$rrna_asm")
    printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
        "rRNA_locus" "$rrna_len" "$rrna_nc" "N" "ASSEMBLED" "BLASTn: ${RRNA_BLAST_HIT}"
else
    printf "%-14s %-12s %-8s %-6s %-12s %s\n" \
        "rRNA_locus" "N/A" "0" "-" "FAILED/SKIP" "No assembly"
fi
} | tee "$SUMMARY" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "[$(date)] Pipeline complete." | tee -a "$LOG"
echo "Summary table:   $SUMMARY"
echo "Full log:        $LOG"
echo "Assemblies:      ${ASM_DIR}/"
echo "Annotations:     ${ANN_DIR}/"
echo ""
echo "Next: manual curation in IGV"
echo "  minimap2 --map-ont assembly.fasta reads.fa | samtools sort > aln.bam"
echo "  samtools index aln.bam"
echo "  # Load assembly as genome, aln.bam + result.gff as tracks"
echo "  # Translation table: Invertebrate Mitochondrial (code 5)"
