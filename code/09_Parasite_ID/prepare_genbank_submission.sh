#!/usr/bin/env bash
# =============================================================================
# prepare_genbank_submission.sh
# Prepare NCBI GenBank submission files for Eustrongylides sp.
# =============================================================================
# Generates submission-ready FASTA, source modifier TSV, and feature table
# for two loci:
#
#   1. ITS region  (ITS1 + 5.8S + ITS2) — rRNA-ITS portal
#   2. COX1        (protein-coding CDS)  — regular nucleotide portal
#
# Also runs quality checks on COX1: frameshifts and internal stop codons
# using invertebrate mitochondrial genetic code 5.
#
# Hardcoded metadata (fixed for this dataset):
#   Organism      : Eustrongylides sp.
#   Host          : Paramormyrops kingsleyae
#   Isolate       : bc1-4_pool
#   Country       : Gabon
#   Collection date: 2019
#
# Usage:
#   bash prepare_genbank_submission.sh [options]
#
# Options:
#   -a  PATH   COX1 assembly FASTA [default: auto-detect]
#   -g  PATH   MITOS2 GFF for COX1  [default: auto-detect]
#   -i  PATH   ITSx output directory [default: auto-detect]
#   -u  STR    ITS unit name to use  [default: unit1]
#   -o  PATH   Output directory      [default: genbank_submission]
#   -h         Show this help
#
# Output layout:
#   $OUTDIR/ITS/
#     Eustrongylides_sp_ITS_region.fa   — FASTA for rRNA-ITS portal
#     source_modifiers.tsv              — source modifiers table
#   $OUTDIR/COX1/
#     Eustrongylides_sp_COX1.fa         — FASTA for nucleotide portal
#     Eustrongylides_sp_COX1.tbl        — feature table with CDS annotation
#     source_modifiers.tsv              — source modifiers table
#   $OUTDIR/check_report.txt            — quality check results
#
# Submission portals:
#   ITS  : https://submit.ncbi.nlm.nih.gov/ → Nucleotide → rRNA/ITS
#            Template: "contains rRNA-ITS region"
#   COX1 : https://submit.ncbi.nlm.nih.gov/ → Nucleotide → Sequences
#            (protein-coding gene — NOT the rRNA portal)
# =============================================================================

set -euo pipefail

# =============================================================================
# Hardcoded metadata
# =============================================================================
ORGANISM="Eustrongylides sp."
HOST="Paramormyrops kingsleyae"
ISOLATE="bc1-4_pool"
COUNTRY="Gabon"
COLLECTION_DATE="2019"
SEQ_ID_ITS="Eustrongylides_sp_bc14pool_ITS"
SEQ_ID_COX1="Eustrongylides_sp_bc14pool_COX1"

# =============================================================================
# Defaults and argument parsing
# =============================================================================

# Auto-detect paths relative to script location
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RUN_DIR="${ROOT}/output_data/09_Parasite_ID/april06_full_run"

COX1_ASM="${RUN_DIR}/assemblies/mtDNA_02/assembly.fasta"
COX1_GFF="${RUN_DIR}/annotations/mtDNA_02/result.gff"
ITSX_DIR="${RUN_DIR}/rrna_annotation/itsx"
ITS_UNIT="unit1"
OUTDIR="genbank_submission"

usage() {
    sed -n '9,27p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 1
}

while getopts "a:g:i:u:o:h" opt; do
    case $opt in
        a) COX1_ASM="$OPTARG" ;;
        g) COX1_GFF="$OPTARG" ;;
        i) ITSX_DIR="$OPTARG" ;;
        u) ITS_UNIT="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

mkdir -p "${OUTDIR}/ITS" "${OUTDIR}/COX1"
REPORT="${OUTDIR}/check_report.txt"
: > "$REPORT"

log() { echo "$@" | tee -a "$REPORT"; }

log "GenBank Submission Preparation"
log "Generated: $(date)"
log "Organism:  ${ORGANISM}"
log "Isolate:   ${ISOLATE}"
log "Country:   ${COUNTRY}  Collection: ${COLLECTION_DATE}"
log "Host:      ${HOST}"
log ""

# =============================================================================
# STEP 1: ITS region — extract, validate, format
# =============================================================================

log "================================================================="
log "STEP 1: ITS region (ITS1 + 5.8S + ITS2)"
log "================================================================="

ITS1_FA="${ITSX_DIR}/${ITS_UNIT}.ITS1.fasta"
S58_FA="${ITSX_DIR}/${ITS_UNIT}.5_8S.fasta"
ITS2_FA="${ITSX_DIR}/${ITS_UNIT}.ITS2.fasta"

for f in "$ITS1_FA" "$S58_FA" "$ITS2_FA"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: ITSx output not found or empty: $f" | tee -a "$REPORT"
        exit 1
    fi
done

ITS_OUT="${OUTDIR}/ITS/Eustrongylides_sp_ITS_region.fa"
ITS_MODS="${OUTDIR}/ITS/source_modifiers.tsv"

python3 - <<PYEOF 2>&1 | tee -a "$REPORT"
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_seq(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return str(rec.seq).upper().replace("-", "")

its1 = read_seq("$ITS1_FA")
s58  = read_seq("$S58_FA")
its2 = read_seq("$ITS2_FA")
merged = its1 + s58 + its2

# Wrap at 60 chars
wrapped = "\n".join(merged[i:i+60] for i in range(0, len(merged), 60))

header = ('>${SEQ_ID_ITS} '
          '[organism=${ORGANISM}] '
          '[isolate=${ISOLATE}] '
          '[country=${COUNTRY}] '
          '[collection_date=${COLLECTION_DATE}] '
          '[host=${HOST}] '
          '[mol_type=genomic DNA]')

with open("$ITS_OUT", "w") as fh:
    fh.write(header + "\n" + wrapped + "\n")

# Quality checks
n_count = merged.count("N") + merged.count("n")
ambig = sum(1 for b in merged if b not in "ACGTacgt")

print(f"ITS1 length:  {len(its1)} bp")
print(f"5.8S length:  {len(s58)} bp")
print(f"ITS2 length:  {len(its2)} bp")
print(f"Total:        {len(merged)} bp")
print(f"N bases:      {n_count}")
print(f"Ambiguous:    {ambig}")
if ambig > 0:
    print("WARNING: ambiguous bases present — NCBI may flag these")
else:
    print("Status: READY")
PYEOF

# Source modifiers TSV
{
printf "Sequence_ID\torganism\tisolate\tcountry\tcollection_date\thost\tmol_type\n"
printf "%s\t%s\t%s\t%s\t%s\t%s\tgenomic DNA\n" \
    "$SEQ_ID_ITS" "$ORGANISM" "$ISOLATE" "$COUNTRY" "$COLLECTION_DATE" "$HOST"
} > "$ITS_MODS"

log ""
log "Output: ${ITS_OUT}"
log "Source modifiers: ${ITS_MODS}"

# =============================================================================
# STEP 2: COX1 CDS — extract, validate, check translation
# =============================================================================

log ""
log "================================================================="
log "STEP 2: COX1 coding sequence"
log "================================================================="

COX1_OUT="${OUTDIR}/COX1/Eustrongylides_sp_COX1.fa"
COX1_TBL="${OUTDIR}/COX1/Eustrongylides_sp_COX1.tbl"
COX1_MODS="${OUTDIR}/COX1/source_modifiers.tsv"

python3 - <<PYEOF 2>&1 | tee -a "$REPORT"
import sys
from Bio import SeqIO
from Bio.Seq import Seq

# ---- Parse MITOS2 GFF: find best cox1 exons ----
# Strategy: collect all exon features for cox1_0 (highest-scoring gene)
# and cox1_1. Use cox1_0 (score ~812M vs ~203M for cox1_1).
# cox1_0 has two exons; join them.

gff_file  = "$COX1_GFF"
asm_file  = "$COX1_ASM"
out_fa    = "$COX1_OUT"
out_tbl   = "$COX1_TBL"
seq_id    = "$SEQ_ID_COX1"
organism  = "$ORGANISM"
isolate   = "$ISOLATE"
country   = "$COUNTRY"
coll_date = "$COLLECTION_DATE"
host      = "$HOST"

exons = {}   # gene_name -> list of (start, end, score) 1-based inclusive
gaps  = []   # gaps between exons of the same gene (for reporting)

with open(gff_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "exon":
            continue
        start  = int(parts[3])
        end    = int(parts[4])
        score  = float(parts[5]) if parts[5] != "." else 0.0
        strand = parts[6]
        attrs  = {a.split("=")[0]: a.split("=")[1]
                  for a in parts[8].split(";") if "=" in a}
        parent = attrs.get("Parent", "")
        name   = attrs.get("Name", "")
        # Group by parent transcript
        if parent not in exons:
            exons[parent] = []
        exons[parent].append({"start": start, "end": end,
                               "score": score, "name": name, "strand": strand})

# Pick the parent with highest total score
best_parent = max(exons, key=lambda p: sum(e["score"] for e in exons[p]))
best_exons  = sorted(exons[best_parent], key=lambda e: e["start"])
strand      = best_exons[0]["strand"]

print(f"Best COX1 gene: {best_parent}")
print(f"Strand: {strand}")
for i, ex in enumerate(best_exons, 1):
    print(f"  Exon {i}: {ex['start']}-{ex['end']} ({ex['end']-ex['start']+1} bp)  score={ex['score']:.1f}")

# Report any gaps between exons
for i in range(1, len(best_exons)):
    gap_start = best_exons[i-1]["end"] + 1
    gap_end   = best_exons[i]["start"] - 1
    gap_len   = gap_end - gap_start + 1
    if gap_len > 0:
        print(f"  Gap between exon {i} and {i+1}: positions {gap_start}-{gap_end} ({gap_len} bp) — excluded from CDS")
        gaps.append((gap_start, gap_end, gap_len))

# ---- Extract and join exon sequences ----
asm_rec = next(SeqIO.parse(asm_file, "fasta"))
asm_seq = str(asm_rec.seq).upper()

joined = ""
for ex in best_exons:
    # GFF is 1-based inclusive → Python 0-based slice
    joined += asm_seq[ex["start"]-1 : ex["end"]]

if strand == "-":
    joined = str(Seq(joined).reverse_complement())

print(f"\nJoined CDS length: {len(joined)} bp")

# ---- Quality checks ----
warnings = []

# Frame check
if len(joined) % 3 != 0:
    msg = f"FRAMESHIFT WARNING: CDS length {len(joined)} is not divisible by 3 (remainder {len(joined)%3})"
    print(msg)
    warnings.append(msg)
else:
    print(f"Frame check: {len(joined)} mod 3 = 0 ✓")

# Translation with invertebrate mitochondrial code (table 5)
# Split into codons and translate manually to catch all stops
codon_table_5_stops = {"TAA", "TAG"}  # TGA = Trp in code 5
codon_table_5_start = {"ATT", "ATC", "ATA", "ATG", "GTG", "TTG", "CTG"}

start_codon = joined[0:3]
stop_codon  = joined[-3:] if len(joined) >= 3 else "???"
aa_len      = len(joined) // 3

print(f"Start codon: {start_codon}", end="")
if start_codon in codon_table_5_start:
    print(" (valid start in code 5 ✓)")
else:
    msg = f"WARNING: unexpected start codon {start_codon}"
    print(f" ← {msg}")
    warnings.append(msg)

# Check for internal stop codons (all triplets except the last)
internal_stops = []
for i in range(0, len(joined)-3, 3):
    codon = joined[i:i+3]
    if codon in codon_table_5_stops:
        aa_pos = i // 3 + 1
        internal_stops.append((aa_pos, i+1, codon))

if internal_stops:
    for aa_pos, nt_pos, codon in internal_stops:
        msg = f"STOP CODON WARNING: {codon} at nt {nt_pos} (aa position {aa_pos}/{aa_len})"
        print(msg)
        warnings.append(msg)
else:
    print(f"Internal stop codons: none ✓")

print(f"Terminal codon: {stop_codon}", end="")
if stop_codon in codon_table_5_stops:
    print(" (stop ✓)")
else:
    print(f" (note: {stop_codon} is not TAA/TAG — may be incomplete CDS)")

print(f"\nCOX1 Status: {'WARNINGS — review before submission' if warnings else 'READY'}")

# ---- Write submission FASTA ----
header = (f">{seq_id} "
          f"[organism={organism}] "
          f"[isolate={isolate}] "
          f"[country={country}] "
          f"[collection_date={coll_date}] "
          f"[host={host}] "
          f"[mol_type=genomic DNA]")
wrapped = "\n".join(joined[i:i+60] for i in range(0, len(joined), 60))
with open(out_fa, "w") as fh:
    fh.write(header + "\n" + wrapped + "\n")

# ---- Write feature table ----
# The submitted FASTA is the already-joined CDS (no intron present in the
# sequence). Use a single <1..>N interval marked partial at both ends.
# A multi-interval feature table is only correct when the submitted sequence
# contains the intron — here it does not, so a split annotation would make
# NCBI look for an intron that isn't there and flag a frameshift.
cds_len = len(joined)
with open(out_tbl, "w") as fh:
    fh.write(f">Feature {seq_id}\n")
    fh.write(f"<1\t>{cds_len}\tCDS\n")
    fh.write(f"\t\t\tproduct\tcytochrome c oxidase subunit 1\n")
    fh.write(f"\t\t\tcodon_start\t1\n")
    fh.write(f"\t\t\ttransl_table\t5\n")
    fh.write(f"\t\t\tprotein_id\tgnl|bc14pool|cox1\n")
    if gaps:
        note = "; ".join(
            f"assembly gap of {g[2]}bp at original contig pos {g[0]}-{g[1]} excluded from joined CDS"
            for g in gaps)
        fh.write(f"\t\t\tnote\t{note}\n")

print(f"\nOutput FASTA:   {out_fa}")
print(f"Feature table:  {out_tbl}")
PYEOF

# Source modifiers TSV
{
printf "Sequence_ID\torganism\tisolate\tcountry\tcollection_date\thost\tmol_type\n"
printf "%s\t%s\t%s\t%s\t%s\t%s\tgenomic DNA\n" \
    "$SEQ_ID_COX1" "$ORGANISM" "$ISOLATE" "$COUNTRY" "$COLLECTION_DATE" "$HOST"
} > "$COX1_MODS"

log ""
log "Output: ${COX1_OUT}"
log "Feature table: ${COX1_TBL}"
log "Source modifiers: ${COX1_MODS}"

# =============================================================================
# STEP 3: Submission summary
# =============================================================================

log ""
log "================================================================="
log "STEP 3: Submission file summary"
log "================================================================="
log ""
log "Submission 1 — rRNA-ITS region"
log "  Portal:   https://submit.ncbi.nlm.nih.gov/ → Nucleotide → rRNA/ITS"
log "  Template: 'contains rRNA-ITS region (ITS1, 5.8S, ITS2)'"
log "  FASTA:    ${OUTDIR}/ITS/Eustrongylides_sp_ITS_region.fa"
log "  Modifiers:${OUTDIR}/ITS/source_modifiers.tsv"
log ""
log "Submission 2 — COX1 (cytochrome c oxidase subunit 1)"
log "  Portal:   https://submit.ncbi.nlm.nih.gov/ → Nucleotide → Sequences"
log "  NOTE:     COX1 is protein-coding — use the regular nucleotide portal,"
log "            NOT the rRNA/ITS portal"
log "  FASTA:    ${OUTDIR}/COX1/Eustrongylides_sp_COX1.fa"
log "  Feature:  ${OUTDIR}/COX1/Eustrongylides_sp_COX1.tbl"
log "  Modifiers:${OUTDIR}/COX1/source_modifiers.tsv"
log ""
log "Full check report: ${REPORT}"
log ""
log "Done."
