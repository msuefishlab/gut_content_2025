#!/usr/bin/env bash
# =============================================================================
# prepare_bankit_submission.sh
# Prepare NCBI BankIt submission files for Eustrongylides sp.
# mitochondrial minichromosomes (mtDNA_01 – mtDNA_10).
# =============================================================================
# Generates three submission-ready files:
#
#   sequences.fasta      — all 10 contigs with NCBI-style inline modifiers
#   feature_table.tbl    — combined NCBI 5-column feature table (all contigs)
#   source_modifiers.tsv — metadata table (one row per contig)
#   check_report.txt     — per-CDS QC (frame, stop codons, ambiguous bases)
#
# Hardcoded metadata (fixed for this dataset):
#   Organism      : Eustrongylides sp.
#   Isolate       : bc1-4_pool
#   Host          : Paramormyrops kingsleyae
#   Country       : Gabon
#   Collection    : 2019
#   Organelle     : mitochondrion
#   Genetic code  : 5  (invertebrate mitochondrial)
#   Topology      : circular
#
# BankIt submission portal:
#   https://www.ncbi.nlm.nih.gov/WebSub/
#   → Genome → Organelle → Circular
#
# Usage:
#   bash prepare_bankit_submission.sh [options]
#
# Options:
#   -a  PATH   Base directory for assemblies   [default: auto-detect]
#   -n  PATH   Base directory for annotations  [default: auto-detect]
#   -o  PATH   Output directory                [default: bankit_submission]
#   -h         Show this help
#
# =============================================================================

set -euo pipefail

# =============================================================================
# Hardcoded metadata
# =============================================================================
ORGANISM="Eustrongylides sp."
ISOLATE="bc1-4_pool"
HOST="Paramormyrops kingsleyae"
COUNTRY="Gabon"
COLLECTION_DATE="2019"
MOL_TYPE="genomic DNA"
ORGANELLE="mitochondrion"
TOPOLOGY="circular"
TRANSL_TABLE=5

# =============================================================================
# Defaults and argument parsing
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RUN_DIR="${ROOT}/output_data/09_Parasite_ID/april06_full_run"

ASM_BASE="${RUN_DIR}/assemblies"
ANN_BASE="${RUN_DIR}/annotations"
OUTDIR="${ROOT}/output_data/09_Parasite_ID/bankit_submission"

usage() {
    sed -n '29,36p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 1
}

while getopts "a:n:o:h" opt; do
    case $opt in
        a) ASM_BASE="$OPTARG" ;;
        n) ANN_BASE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

mkdir -p "${OUTDIR}"
REPORT="${OUTDIR}/check_report.txt"
: > "$REPORT"

log() { echo "$@" | tee -a "$REPORT"; }

log "BankIt Submission Preparation — Mitochondrial Minichromosomes"
log "Generated:    $(date)"
log "Organism:     ${ORGANISM}"
log "Isolate:      ${ISOLATE}"
log "Host:         ${HOST}"
log "Country:      ${COUNTRY}  Collection: ${COLLECTION_DATE}"
log "Organelle:    ${ORGANELLE}  Topology: ${TOPOLOGY}"
log "Genetic code: ${TRANSL_TABLE} (invertebrate mitochondrial)"
log ""

# =============================================================================
# STEP 1: Build combined FASTA with NCBI inline modifiers
# =============================================================================

log "================================================================="
log "STEP 1: Building combined FASTA"
log "================================================================="

FASTA_OUT="${OUTDIR}/sequences.fasta"
: > "$FASTA_OUT"

# Per-chromosome leading trim (0-based start index).
# Trims the specified number of bases from the start of the assembly before
# submission. For circular sequences this just shifts the origin.
# GFF coordinate offsets are applied in Step 2 to match.
chr_trim() {
    # Usage: chr_trim CHR_ID → prints trim offset (0 if not trimmed)
    case "$1" in
        mtDNA_07) echo 112 ;;   # remove 0-111 (VecScreen moderate hit)
        *)        echo 0   ;;
    esac
}

for CHR_NUM in 01 02 03 04 05 06 07 08 09 10; do
    CHR="mtDNA_${CHR_NUM}"
    ASM="${ASM_BASE}/${CHR}/assembly.fasta"

    if [[ ! -s "$ASM" ]]; then
        log "  SKIP ${CHR}: assembly not found at ${ASM}"
        continue
    fi

    SEQ_ID="Eustrongylides_sp_bc14pool_${CHR}"
    TRIM="$(chr_trim "${CHR}")"

    python3 - <<PYEOF >> "$FASTA_OUT"
from Bio import SeqIO

# When assembly has multiple contigs, use the longest one (the true
# minichromosome; shorter contigs are assembly artifacts).
recs = list(SeqIO.parse("${ASM}", "fasta"))
rec = max(recs, key=lambda r: len(r.seq))
seq = str(rec.seq).upper()
trim = int("${TRIM}")
if trim:
    seq = seq[trim:]

header = (
    ">${SEQ_ID} "
    "[organism=${ORGANISM}] "
    "[isolate=${ISOLATE}] "
    "[host=${HOST}] "
    "[country=${COUNTRY}] "
    "[collection_date=${COLLECTION_DATE}] "
    "[mol_type=${MOL_TYPE}]"
    # organelle and topology are set via source_modifiers.tsv and the BankIt form,
    # NOT as inline FASTA modifiers (BankIt does not recognise them in headers)
)
wrapped = "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))
print(header)
print(wrapped)
if trim:
    import sys
    print(f"  {chr_id if False else '${CHR}'}: trimmed first {trim} bases (new length {len(seq)} bp)", file=sys.stderr)
PYEOF

    log "  ${CHR}: trimmed_start=${TRIM} → ${SEQ_ID}"
done

log ""
log "Combined FASTA: ${FASTA_OUT}"
log ""

# =============================================================================
# STEP 2 + 3: Convert MITOS2 GFF → feature table; run CDS quality checks
# =============================================================================

log "================================================================="
log "STEP 2: Building feature table + CDS quality checks"
log "================================================================="

TBL_OUT="${OUTDIR}/feature_table.tbl"
: > "$TBL_OUT"

# Export variables needed by the embedded Python script
export ASM_BASE ANN_BASE TBL_OUT ISOLATE REPORT

python3 - <<'PYEOF' 2>&1 | tee -a "$REPORT"
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Configuration (injected at runtime by the bash heredoc)
# ---------------------------------------------------------------------------
ASM_BASE   = os.environ.get("ASM_BASE",  "")
ANN_BASE   = os.environ.get("ANN_BASE",  "")
TBL_OUT    = os.environ.get("TBL_OUT",   "")
ISOLATE    = os.environ.get("ISOLATE",   "bc1-4_pool")
REPORT_FH  = open(os.environ.get("REPORT", "/dev/null"), "a")

CHR_NUMS   = [f"{n:02d}" for n in range(1, 11)]

# Per-chromosome leading trim (0-based bases removed from start of assembly).
# Must match CHR_TRIM_START in the bash section above.
# All GFF coordinates are shifted left by this amount.
CHR_TRIM = {
    "mtDNA_07": 112,   # remove 0-111 (VecScreen moderate hit)
}

# ---------------------------------------------------------------------------
# Name mappings: MITOS2 gene name → (ncbi_gene_tag, ncbi_product)
# ---------------------------------------------------------------------------
GENE_NAMES = {
    "atp6":  ("atp6",  "ATP synthase F0 subunit 6"),
    "cox1":  ("cox1",  "cytochrome c oxidase subunit 1"),
    "cox2":  ("cox2",  "cytochrome c oxidase subunit 2"),
    "cox3":  ("cox3",  "cytochrome c oxidase subunit 3"),
    "cytb":  ("cytb",  "cytochrome b"),
    "cob":   ("cytb",  "cytochrome b"),   # MITOS2 alternate name for cytb
    "nad1":  ("nad1",  "NADH dehydrogenase subunit 1"),
    "nad2":  ("nad2",  "NADH dehydrogenase subunit 2"),
    "nad3":  ("nad3",  "NADH dehydrogenase subunit 3"),
    "nad4":  ("nad4",  "NADH dehydrogenase subunit 4"),
    "nad4l": ("nad4l", "NADH dehydrogenase subunit 4L"),
    "nad5":  ("nad5",  "NADH dehydrogenase subunit 5"),
    "nad6":  ("nad6",  "NADH dehydrogenase subunit 6"),
}

TRNA_AA = {
    "trna":  "Ala",  # fallback
    "trnA":  "Ala", "trnC": "Cys", "trnD": "Asp", "trnE": "Glu",
    "trnF":  "Phe", "trnG": "Gly", "trnH": "His", "trnI": "Ile",
    "trnK":  "Lys", "trnL1": "Leu", "trnL2": "Leu",
    "trnM":  "Met", "trnN": "Asn", "trnP": "Pro", "trnQ": "Gln",
    "trnR":  "Arg", "trnS1": "Ser", "trnS2": "Ser",
    "trnT":  "Thr", "trnV": "Val", "trnW": "Trp", "trnY": "Tyr",
}

RRNA_PRODUCTS = {
    "rrn12": "12S ribosomal RNA",
    "rrn16": "16S ribosomal RNA",
    "rrn5":  "5S ribosomal RNA",
}

# Code-5 rules for CDS validation
CODE5_STOPS = {"TAA", "TAG"}
CODE5_START = {"ATT", "ATC", "ATA", "ATG", "GTG", "TTG", "CTG"}

def parse_attrs(attr_str):
    """Parse GFF3 attribute string into a dict."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs

def gff_feature_type(mitos_type, gene_name):
    """Classify a MITOS2 gene-level feature into pcg/trna/rrna/rep_origin."""
    if mitos_type in ("gene",):
        base = gene_name.lower().split("_")[0]
        if base in GENE_NAMES:
            return "pcg"
    if mitos_type == "ncRNA_gene":
        base = gene_name.lower().split("_")[0]
        if base.startswith("trn"):
            return "trna"
        if base.startswith("rrn"):
            return "rrna"
    if mitos_type == "origin_of_replication":
        return "rep_origin"
    return None

def validate_cds(seq_str, gene_name, chr_id, out):
    """Run code-5 CDS quality checks; write results to out."""
    out.append(f"  CDS check for {gene_name} on {chr_id}:")
    warnings = []

    if len(seq_str) % 3 != 0:
        msg = f"    FRAMESHIFT: length {len(seq_str)} not divisible by 3"
        out.append(msg)
        warnings.append(msg)
    else:
        out.append(f"    Frame:        {len(seq_str)} bp (ok)")

    if len(seq_str) < 3:
        out.append("    ERROR: sequence too short to check codons")
        return warnings

    start = seq_str[:3]
    stop  = seq_str[-3:]

    if start in CODE5_START:
        out.append(f"    Start codon:  {start} (valid in code 5)")
    else:
        msg = f"    WARNING: unexpected start codon {start}"
        out.append(msg)
        warnings.append(msg)

    internal_stops = []
    for i in range(0, len(seq_str) - 3, 3):
        codon = seq_str[i:i+3]
        if codon in CODE5_STOPS:
            internal_stops.append((i // 3 + 1, i + 1, codon))

    if internal_stops:
        for aa_pos, nt_pos, codon in internal_stops:
            msg = f"    STOP CODON: {codon} at nt {nt_pos} (aa {aa_pos})"
            out.append(msg)
            warnings.append(msg)
    else:
        out.append("    Internal stops: none (ok)")

    if stop in CODE5_STOPS:
        out.append(f"    Terminal codon: {stop} (stop)")
    else:
        out.append(f"    Terminal codon: {stop} (incomplete CDS or non-standard stop)")

    n_ambig = sum(1 for b in seq_str if b not in "ACGT")
    out.append(f"    Ambiguous bases: {n_ambig}")
    if n_ambig > 0:
        warnings.append(f"    WARNING: {n_ambig} ambiguous bases")

    status = "WARNINGS — review before submission" if warnings else "READY"
    out.append(f"    Status: {status}")
    return warnings

# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
all_warnings = {}

with open(TBL_OUT, "w") as tbl:
    for num in CHR_NUMS:
        chr_id = f"mtDNA_{num}"
        seq_id = f"Eustrongylides_sp_bc14pool_{chr_id}"
        asm_fa = os.path.join(ASM_BASE, chr_id, "assembly.fasta")
        gff_f  = os.path.join(ANN_BASE, chr_id, "result.gff")

        log_lines = [f"\n--- {chr_id} ---"]

        if not os.path.isfile(asm_fa):
            log_lines.append(f"  SKIP: assembly not found at {asm_fa}")
            print("\n".join(log_lines))
            continue

        # Resolve GFF path: try direct first, then numbered subdirs.
        # When subdirs exist, pick the one whose result.gff has the most
        # non-header feature lines (most annotated contig).
        ann_dir = os.path.join(ANN_BASE, chr_id)
        if os.path.isfile(gff_f):
            best_gff = gff_f
            gff_contig = None   # single-contig assembly; no filtering needed
        else:
            candidates = []
            for entry in sorted(os.listdir(ann_dir)):
                candidate = os.path.join(ann_dir, entry, "result.gff")
                if os.path.isfile(candidate):
                    with open(candidate) as fh:
                        n_features = sum(
                            1 for l in fh
                            if not l.startswith("#") and len(l.strip().split("\t")) >= 9
                        )
                    candidates.append((n_features, entry, candidate))
            if not candidates:
                log_lines.append(f"  SKIP: no result.gff found under {ann_dir}")
                print("\n".join(log_lines))
                continue
            candidates.sort(reverse=True)
            best_n, best_subdir, best_gff = candidates[0]
            log_lines.append(f"  NOTE: using subdir '{best_subdir}' ({best_n} features); skipped: {[c[1] for c in candidates[1:]]}")
            # Determine which contig this subdir annotated (from region feature)
            gff_contig = None
            with open(best_gff) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3 and parts[2] == "region":
                        gff_contig = parts[0]
                        break

        # Load assembly sequence; select the correct contig if specified
        trim_offset = CHR_TRIM.get(chr_id, 0)
        asm_recs = {r.id: r for r in SeqIO.parse(asm_fa, "fasta")}
        if gff_contig and gff_contig in asm_recs:
            asm_seq = str(asm_recs[gff_contig].seq).upper()[trim_offset:]
            log_lines.append(f"  NOTE: using contig '{gff_contig}' from assembly")
        else:
            # Default: first contig
            asm_seq = str(next(iter(asm_recs.values())).seq).upper()[trim_offset:]
        if trim_offset:
            log_lines.append(f"  NOTE: applying trim offset -{trim_offset} to all coordinates")

        gff_f = best_gff  # use resolved path for parsing below

        # ----------------------------------------------------------------
        # Parse GFF: collect gene-level features and exons
        # ----------------------------------------------------------------
        gene_features  = []   # (start, end, strand, mitos_type, gene_name, gene_id)
        exon_by_parent = {}   # parent_transcript → [(start, end, score, strand)]
        rep_origins    = []   # (start, end, strand, name)

        with open(gff_f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                feat_type = parts[2]
                start  = int(parts[3]) - trim_offset
                end    = int(parts[4]) - trim_offset
                # Skip features that fall entirely within the trimmed region
                if end <= 0:
                    continue
                score  = float(parts[5]) if parts[5] not in (".", "") else 0.0
                strand = parts[6]
                attrs  = parse_attrs(parts[8])
                name   = attrs.get("Name", "")
                gene_id = attrs.get("gene_id", name)

                if feat_type in ("gene", "ncRNA_gene"):
                    gene_features.append((start, end, strand, feat_type, name, gene_id))

                elif feat_type == "exon":
                    parent = attrs.get("Parent", "")
                    if parent not in exon_by_parent:
                        exon_by_parent[parent] = []
                    exon_by_parent[parent].append(
                        {"start": start, "end": end, "score": score, "strand": strand}
                    )

                elif feat_type == "origin_of_replication":
                    rep_origins.append((start, end, strand, name))

        # ----------------------------------------------------------------
        # Write feature table block
        # ----------------------------------------------------------------
        tbl.write(f">Feature {seq_id}\n")
        chr_warnings = []

        # ---- Protein-coding and RNA gene features ----
        # Build parent → total score mapping
        parent_scores = {}
        for parent, exons in exon_by_parent.items():
            parent_scores[parent] = sum(e["score"] for e in exons)

        # Deduplicate PCG gene features: MITOS2 sometimes reports multiple
        # copies (e.g. cox1_0, cox1_1). Keep only the one whose matching
        # transcript parent has the highest total exon score.
        seen_pcg = {}   # base_name → (g_start, g_end, g_strand, g_type, g_name, g_id, best_score)
        filtered_gene_features = []
        for feat in gene_features:
            g_start, g_end, g_strand, g_type, g_name, g_id = feat
            ftype = gff_feature_type(g_type, g_name)
            if ftype == "pcg":
                base_name = g_name.lower().split("_")[0]
                # Score = highest matching transcript score
                candidate_parents = [
                    p for p in exon_by_parent
                    if p.lower().startswith(f"transcript_{base_name}")
                ]
                best_score = max(
                    (parent_scores.get(p, 0) for p in candidate_parents),
                    default=0
                )
                if base_name not in seen_pcg or best_score > seen_pcg[base_name][1]:
                    seen_pcg[base_name] = (feat, best_score)
            else:
                filtered_gene_features.append(feat)
        # Add deduplicated PCG features (one per gene)
        for base_name, (feat, _score) in seen_pcg.items():
            filtered_gene_features.append(feat)

        # For each gene feature, find matching exon parents
        for (g_start, g_end, g_strand, g_type, g_name, g_id) in filtered_gene_features:
            ftype = gff_feature_type(g_type, g_name)
            if ftype is None:
                continue

            base_name = g_name.lower().split("_")[0]  # e.g. "cox1" from "cox1_0"
            # Canonical name normalisation
            canon_name = g_name.split("_")[0]  # preserve case for tRNA

            # ---- Protein-coding gene (PCG) ----
            if ftype == "pcg":
                gene_tag, product = GENE_NAMES.get(base_name, (base_name, base_name))

                # Find exon parents for this gene (match by gene name prefix)
                candidate_parents = [
                    p for p in exon_by_parent
                    if p.lower().startswith(f"transcript_{base_name}")
                ]
                if not candidate_parents:
                    # Fall back: use gene coordinates as single exon
                    log_lines.append(f"  NOTE: no exon features for {g_name}; using gene bounds")
                    best_exons = [{"start": g_start, "end": g_end, "score": 0.0, "strand": g_strand}]
                else:
                    # Pick parent with highest total score
                    best_parent = max(candidate_parents, key=lambda p: parent_scores.get(p, 0))
                    best_exons  = sorted(exon_by_parent[best_parent], key=lambda e: e["start"])
                    log_lines.append(f"  PCG {gene_tag}: using {best_parent} ({len(best_exons)} exon(s))")

                strand = best_exons[0]["strand"]

                # Run QC first so we know whether CDS is safe to annotate
                joined = ""
                for s, e in [(ex["start"], ex["end"]) for ex in best_exons]:
                    joined += asm_seq[s-1:e]
                if strand == "-":
                    joined = str(Seq(joined).reverse_complement())
                w = validate_cds(joined, gene_tag, chr_id, log_lines)
                chr_warnings.extend(w)
                has_internal_stops = any(
                    "STOP CODON" in msg or "FRAMESHIFT" in msg for msg in w
                )

                # Build coordinate lists (NCBI: first row has feature key,
                # subsequent rows blank; minus-strand listed high→low)
                e_coords = [(e["start"], e["end"]) for e in best_exons]
                if strand == "-":
                    e_coords_r = list(reversed([(e, s) for s, e in e_coords]))
                else:
                    e_coords_r = None   # not used for plus strand

                # Write gene feature
                coord_list = e_coords_r if strand == "-" else e_coords
                for i, (s, e) in enumerate(coord_list):
                    tbl.write(f"{s}\t{e}\t{'gene' if i == 0 else ''}\n")
                tbl.write(f"\t\t\tgene\t{gene_tag}\n")

                if has_internal_stops:
                    # CDS has internal stops — omit CDS feature to avoid
                    # BankIt rejection; record as a note on the gene instead.
                    log_lines.append(f"  WARNING: CDS for {gene_tag} suppressed (internal stop codons detected; MITOS2 annotation may be incorrect)")
                    chr_warnings.append(f"CDS suppressed for {gene_tag} on {chr_id} due to internal stop codons")
                else:
                    # Write CDS feature
                    for i, (s, e) in enumerate(coord_list):
                        tbl.write(f"{s}\t{e}\t{'CDS' if i == 0 else ''}\n")
                    tbl.write(f"\t\t\tproduct\t{product}\n")
                    tbl.write(f"\t\t\tcodon_start\t1\n")
                    tbl.write(f"\t\t\ttransl_table\t5\n")
                    tbl.write(f"\t\t\tprotein_id\tgnl|bc14pool|{gene_tag}_{chr_id}\n")

            # ---- tRNA ----
            elif ftype == "trna":
                # Look up three-letter amino acid
                aa = TRNA_AA.get(canon_name, TRNA_AA.get(base_name, "???"))
                log_lines.append(f"  tRNA {canon_name}: tRNA-{aa}")

                if g_strand == "-":
                    tbl.write(f"{g_end}\t{g_start}\tgene\n")
                else:
                    tbl.write(f"{g_start}\t{g_end}\tgene\n")
                tbl.write(f"\t\t\tgene\t{canon_name}\n")

                if g_strand == "-":
                    tbl.write(f"{g_end}\t{g_start}\ttRNA\n")
                else:
                    tbl.write(f"{g_start}\t{g_end}\ttRNA\n")
                tbl.write(f"\t\t\tproduct\ttRNA-{aa}\n")

            # ---- rRNA ----
            elif ftype == "rrna":
                product = RRNA_PRODUCTS.get(base_name, base_name)
                log_lines.append(f"  rRNA {canon_name}: {product}")

                if g_strand == "-":
                    tbl.write(f"{g_end}\t{g_start}\tgene\n")
                else:
                    tbl.write(f"{g_start}\t{g_end}\tgene\n")
                tbl.write(f"\t\t\tgene\t{canon_name}\n")

                if g_strand == "-":
                    tbl.write(f"{g_end}\t{g_start}\trRNA\n")
                else:
                    tbl.write(f"{g_start}\t{g_end}\trRNA\n")
                tbl.write(f"\t\t\tproduct\t{product}\n")

        # ---- rep_origin features ----
        for (ro_start, ro_end, ro_strand, ro_name) in rep_origins:
            log_lines.append(f"  rep_origin: {ro_name} {ro_start}-{ro_end} ({ro_strand})")
            if ro_strand == "-":
                tbl.write(f"{ro_end}\t{ro_start}\trep_origin\n")
            else:
                tbl.write(f"{ro_start}\t{ro_end}\trep_origin\n")
            tbl.write(f"\t\t\tnote\tputative mitochondrial control region / origin of replication\n")

        all_warnings[chr_id] = chr_warnings
        print("\n".join(log_lines))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n\n=================================================================")
print("QC SUMMARY")
print("=================================================================")
any_warn = False
for chr_id, warns in all_warnings.items():
    if warns:
        any_warn = True
        print(f"  {chr_id}: {len(warns)} warning(s)")
        for w in warns:
            print(f"    {w.strip()}")
    else:
        print(f"  {chr_id}: OK")

if not any_warn:
    print("\nAll CDS passed quality checks — ready for submission.")
else:
    print("\nReview warnings above before submitting to BankIt.")
PYEOF

log ""
log "Feature table: ${TBL_OUT}"
log ""

# =============================================================================
# STEP 4: Source modifier table
# =============================================================================

log "================================================================="
log "STEP 4: Source modifier table"
log "================================================================="

SRC_OUT="${OUTDIR}/source_modifiers.tsv"

{
printf "Sequence_ID\torganism\tisolate\thost\tcountry\tcollection_date\tmol_type\torganelle\ttopology\n"
for CHR_NUM in 01 02 03 04 05 06 07 08 09 10; do
    CHR="mtDNA_${CHR_NUM}"
    SEQ_ID="Eustrongylides_sp_bc14pool_${CHR}"
    ASM="${ASM_BASE}/${CHR}/assembly.fasta"
    [[ -s "$ASM" ]] || continue
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SEQ_ID" "$ORGANISM" "$ISOLATE" "$HOST" \
        "$COUNTRY" "$COLLECTION_DATE" "$MOL_TYPE" "$ORGANELLE" "$TOPOLOGY"
done
} > "$SRC_OUT"

log "Source modifiers: ${SRC_OUT}"
log ""

# =============================================================================
# STEP 5: Submission summary
# =============================================================================

log "================================================================="
log "STEP 5: BankIt submission summary"
log "================================================================="
log ""
log "Files to upload to BankIt:"
log "  FASTA:          ${FASTA_OUT}"
log "  Feature table:  ${TBL_OUT}"
log "  Source mods:    ${SRC_OUT}"
log "  QC report:      ${REPORT}"
log ""
log "BankIt portal:    https://www.ncbi.nlm.nih.gov/WebSub/"
log "  → Genome → Organelle → Circular"
log "  → Genetic code: Invertebrate Mitochondrial (table 5)"
log ""
log "Quick validation (requires table2asn on HPCC):"
log "  table2asn -i sequences.fasta -f feature_table.tbl \\"
log "            -src source_modifiers.tsv -V vb"
log ""
log "Done."
