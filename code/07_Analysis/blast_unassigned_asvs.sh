#!/bin/bash
#
# Script to extract high-abundance unassigned ASVs and prepare them for BLAST against NCBI nt
#
# This script:
# 1. Exports ASV sequences from QIIME2 artifacts
# 2. Extracts high-abundance unassigned ASVs from each primer set
# 3. Creates FASTA files for BLAST analysis
#
# Usage: bash blast_unassigned_asvs.sh [min_reads] [n_samples]
#   min_reads: Minimum read count to consider (default: 100)
#   n_samples: Number of ASVs to sample per primer set (default: 75)

root="$(git rev-parse --show-toplevel)"
source ${root}/gut_contents.env

# Parameters
MIN_READS=${1:-100}  # Default minimum 100 reads (focus on high-abundance)
N_SAMPLES=${2:-75}   # Default sample 75 ASVs per primer set

# Output directories
BLAST_DIR="${root}/output_data/07_Analysis/blast_unassigned"
SEQ_DIR="${BLAST_DIR}/sequences"
EXPORT_DIR="${BLAST_DIR}/exported_seqs"
TMPDIR="$SCRATCH/qiime2_tmp"

mkdir -p ${TMPDIR}

mkdir -p ${SEQ_DIR}
mkdir -p ${EXPORT_DIR}

echo "=================================================="
echo "Extracting unassigned ASVs for BLAST analysis"
echo "=================================================="
echo "Minimum reads threshold: ${MIN_READS}"
echo "Samples per primer set: ${N_SAMPLES}"
echo ""

# Export sequences from QIIME2 artifacts if not already done
for PRIMERSET in primerset1 primerset2; do
    echo "Processing ${PRIMERSET}..."

    QZA_FILE="${root}/output_data/03_Clustered_Data/${PRIMERSET}_all_p985_seqs.qza"
    EXPORT_SUBDIR="${EXPORT_DIR}/${PRIMERSET}"

    if [ ! -f "${EXPORT_SUBDIR}/dna-sequences.fasta" ]; then
        echo "  Exporting sequences from QIIME2 artifact..."
        singularity exec --bind $TMPDIR:/home/qiime2/q2cli \
            ${qiime_image} \
            qiime tools export \
            --input-path ${QZA_FILE} \
            --output-path ${EXPORT_SUBDIR}
    else
        echo "  Sequences already exported, skipping..."
    fi
done

echo ""
echo "Creating subset FASTA files for high-abundance unassigned ASVs..."

# Run R script to extract and sample sequences
singularity exec ${rimage} Rscript - <<'RSCRIPT'
library(tidyverse)
library(Biostrings)

# Get command line args from environment
min_reads <- as.numeric(Sys.getenv("MIN_READS"))
n_samples <- as.numeric(Sys.getenv("N_SAMPLES"))
root <- Sys.getenv("root")

print($root)

blast_dir <- file.path(root, "output_data/07_Analysis/blast_unassigned")
seq_dir <- file.path(blast_dir, "sequences")
export_dir <- file.path(blast_dir, "exported_seqs")

# Read the removed ASVs list
removed_asvs <- read_csv(
    file.path(root, "output_data/07_Analysis/filtering_loss_analysis/tables/removed_asvs_full_list.csv"),
    show_col_types = FALSE
)

cat("\n")
cat("Summary of unassigned ASVs by primer set:\n")
cat("==========================================\n")
summary_table <- removed_asvs %>%
    filter(Category == "Unassigned") %>%
    group_by(PrimerSet) %>%
    summarise(
        Total_ASVs = n(),
        Total_Reads = sum(Total_Reads),
        Mean_Reads = mean(Total_Reads),
        Median_Reads = median(Total_Reads),
        Max_Reads = max(Total_Reads),
        High_Abundance = sum(Total_Reads >= min_reads)
    )
print(summary_table)
cat("\n")

# Process each primer set
for (ps in c("PS1", "PS2")) {
    cat(sprintf("Processing %s...\n", ps))

    # Filter to high-abundance unassigned ASVs for this primer set
    high_abund <- removed_asvs %>%
        filter(Category == "Unassigned",
               PrimerSet == ps,
               Total_Reads >= min_reads) %>%
        arrange(desc(Total_Reads))

    cat(sprintf("  Found %d unassigned ASVs with >= %d reads\n",
                nrow(high_abund), min_reads))

    if (nrow(high_abund) == 0) {
        cat(sprintf("  WARNING: No ASVs found with >= %d reads, skipping %s\n",
                    min_reads, ps))
        next
    }

    # Sample ASVs (weighted by abundance to prioritize high-abundance ones)
    n_to_sample <- min(n_samples, nrow(high_abund))

    set.seed(42)  # For reproducibility
    sampled_indices <- sample(1:nrow(high_abund),
                             size = n_to_sample,
                             replace = FALSE,
                             prob = high_abund$Total_Reads / sum(high_abund$Total_Reads))

    sampled_asvs <- high_abund[sampled_indices, ]

    cat(sprintf("  Sampled %d ASVs (abundance range: %d - %d reads)\n",
                n_to_sample,
                min(sampled_asvs$Total_Reads),
                max(sampled_asvs$Total_Reads)))

    # Read the corresponding FASTA file
    primerset_name <- ifelse(ps == "PS1", "primerset1", "primerset2")
    fasta_file <- file.path(export_dir, primerset_name, "dna-sequences.fasta")

    if (!file.exists(fasta_file)) {
        cat(sprintf("  ERROR: FASTA file not found: %s\n", fasta_file))
        next
    }

    all_seqs <- readDNAStringSet(fasta_file)

    # Extract sequences for sampled ASVs
    asv_ids <- sampled_asvs$ASV
    matching_seqs <- all_seqs[names(all_seqs) %in% asv_ids]

    cat(sprintf("  Matched %d sequences\n", length(matching_seqs)))

    # Add abundance information to sequence names for reference
    names(matching_seqs) <- sapply(names(matching_seqs), function(asv_id) {
        reads <- sampled_asvs$Total_Reads[sampled_asvs$ASV == asv_id]
        sprintf("%s|reads=%d", asv_id, reads)
    })

    # Write output FASTA
    output_fasta <- file.path(seq_dir, sprintf("%s_unassigned_highAbund.fasta", tolower(ps)))
    writeXStringSet(matching_seqs, output_fasta)

    cat(sprintf("  Wrote %s\n", output_fasta))

    # Save metadata
    metadata_file <- file.path(seq_dir, sprintf("%s_unassigned_metadata.tsv", tolower(ps)))
    write_tsv(sampled_asvs, metadata_file)
    cat(sprintf("  Wrote %s\n\n", metadata_file))
}

cat("Done!\n")
RSCRIPT

echo ""
echo "=================================================="
echo "Preparation complete!"
echo "=================================================="
echo "Output files are in: ${BLAST_DIR}"
echo ""
echo "Next steps:"
echo "1. Review the FASTA files in ${SEQ_DIR}"
echo "2. Submit BLAST jobs using blast_unassigned_submit.sb"
echo ""
