#!/usr/bin/env Rscript
#
# Analyze BLAST results for unassigned ASVs to characterize what's being lost
#
# This script:
# 1. Reads BLAST results for both primer sets
# 2. Categorizes hits (metazoan, non-metazoan, no hit)
# 3. Analyzes taxonomic composition of unassigned ASVs
# 4. Generates summary statistics and visualizations

library(tidyverse)
library(patchwork)

# Set up paths
args <- commandArgs(trailingOnly = TRUE)
root <- ifelse(length(args) > 0, args[1], Sys.getenv("root"))

blast_dir <- file.path(root, "output_data/07_Analysis/blast_unassigned")
results_dir <- file.path(blast_dir, "results")
plots_dir <- file.path(blast_dir, "plots")

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("================================================\n")
cat("Analyzing BLAST Results for Unassigned ASVs\n")
cat("================================================\n\n")

# Read BLAST results for both primer sets
results_list <- list()

for (ps in c("ps1", "ps2")) {
    summary_file <- file.path(results_dir, sprintf("%s_blast_summary.tsv", ps))
    metadata_file <- file.path(blast_dir, "sequences", sprintf("%s_unassigned_metadata.tsv", ps))

    if (!file.exists(summary_file)) {
        cat(sprintf("WARNING: BLAST results not found for %s: %s\n", ps, summary_file))
        next
    }

    # Read BLAST summary
    blast_res <- read_tsv(summary_file, show_col_types = FALSE)

    # Read metadata
    if (file.exists(metadata_file)) {
        metadata <- read_tsv(metadata_file, show_col_types = FALSE)
    } else {
        cat(sprintf("WARNING: Metadata not found for %s\n", ps))
        metadata <- NULL
    }

    # Add primer set column
    blast_res$PrimerSet <- toupper(ps)

    results_list[[ps]] <- list(blast = blast_res, metadata = metadata)
}

# Combine results
all_blast <- bind_rows(lapply(results_list, function(x) x$blast))

if (nrow(all_blast) == 0) {
    cat("ERROR: No BLAST results found. Please run BLAST jobs first.\n")
    quit(status = 1)
}

cat(sprintf("Loaded BLAST results for %d ASVs\n\n", nrow(all_blast)))

# Categorize hits
all_blast <- all_blast %>%
    mutate(
        # Extract major taxonomic groups
        Top_Hit_Domain = case_when(
            str_detect(Top_Hit_Species, "Bacteria|Archaea") ~ "Bacteria/Archaea",
            str_detect(Top_Hit_Species, "Fungi") ~ "Fungi",
            str_detect(Top_Hit_Species, "Viridiplantae") ~ "Plants",
            TRUE ~ "Metazoa/Other"
        ),

        # Categorize quality of hit
        Hit_Quality = case_when(
            Percent_Identity >= 97 ~ "Excellent (≥97%)",
            Percent_Identity >= 90 ~ "Good (90-97%)",
            Percent_Identity >= 80 ~ "Moderate (80-90%)",
            TRUE ~ "Poor (<80%)"
        ),

        # Check if it's likely metazoan (animal)
        Likely_Metazoan = str_detect(Top_Hit_Species,
            "Arthropod|Insect|Crustacea|Mollusca|Annelida|Nematoda|Platyhelminthes|Fish|Amphibia|Reptilia|Aves|Mammalia|Cnidaria|Echinodermata"
        ),

        # Check if it's likely COI
        Likely_COI = str_detect(Top_Hit_Title,
            "cytochrome c oxidase|cytochrome oxidase|COI|CO1|COX1"
        )
    )

# Summary statistics
cat("Summary Statistics by Primer Set\n")
cat("=================================\n\n")

summary_stats <- all_blast %>%
    group_by(PrimerSet) %>%
    summarise(
        N_ASVs = n(),
        Total_Reads = sum(Reads),
        Mean_Identity = mean(Percent_Identity),
        Median_Identity = median(Percent_Identity),
        Excellent_Hits = sum(Hit_Quality == "Excellent (≥97%)"),
        Good_Hits = sum(Hit_Quality == "Good (90-97%)"),
        Moderate_Hits = sum(Hit_Quality == "Moderate (80-90%)"),
        Poor_Hits = sum(Hit_Quality == "Poor (<80%)"),
        Likely_Metazoan = sum(Likely_Metazoan),
        Likely_COI = sum(Likely_COI)
    )

print(summary_stats)
cat("\n")

# Taxonomic breakdown
cat("Taxonomic Breakdown\n")
cat("===================\n\n")

tax_breakdown <- all_blast %>%
    group_by(PrimerSet, Top_Hit_Domain) %>%
    summarise(
        N_ASVs = n(),
        Total_Reads = sum(Reads),
        Pct_Reads = Total_Reads / sum(all_blast$Reads[all_blast$PrimerSet == first(PrimerSet)]) * 100,
        .groups = "drop"
    ) %>%
    arrange(PrimerSet, desc(Total_Reads))

print(tax_breakdown)
cat("\n")

# Top species hits
cat("Top 20 Species Hits (by read abundance)\n")
cat("=======================================\n\n")

top_species <- all_blast %>%
    group_by(PrimerSet, Top_Hit_Species) %>%
    summarise(
        N_ASVs = n(),
        Total_Reads = sum(Reads),
        Mean_Identity = mean(Percent_Identity),
        .groups = "drop"
    ) %>%
    group_by(PrimerSet) %>%
    slice_max(order_by = Total_Reads, n = 20) %>%
    arrange(PrimerSet, desc(Total_Reads))

print(top_species)
cat("\n")

# Save detailed results
output_file <- file.path(results_dir, "blast_analysis_summary.csv")
write_csv(all_blast, output_file)
cat(sprintf("Saved detailed results to: %s\n\n", output_file))

# Create visualizations
cat("Generating plots...\n")

# Plot 1: Hit quality distribution
p1 <- ggplot(all_blast, aes(x = Hit_Quality, y = Reads, fill = PrimerSet)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_log10(labels = scales::comma) +
    labs(
        title = "Unassigned ASV Hit Quality Distribution",
        subtitle = "Total reads per hit quality category",
        x = "Hit Quality (% Identity)",
        y = "Total Reads (log scale)",
        fill = "Primer Set"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Taxonomic composition
p2 <- ggplot(tax_breakdown, aes(x = PrimerSet, y = Total_Reads, fill = Top_Hit_Domain)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = scales::comma) +
    labs(
        title = "Taxonomic Composition of Unassigned ASVs",
        subtitle = "By read abundance",
        x = "Primer Set",
        y = "Total Reads",
        fill = "Taxonomic Group"
    ) +
    theme_bw()

# Plot 3: Identity distribution
p3 <- ggplot(all_blast, aes(x = Percent_Identity, fill = PrimerSet)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
    geom_vline(xintercept = c(80, 90, 97), linetype = "dashed", color = "gray40") +
    facet_wrap(~PrimerSet, ncol = 1) +
    labs(
        title = "BLAST Identity Distribution",
        subtitle = "Dashed lines at 80%, 90%, and 97% identity",
        x = "Percent Identity",
        y = "Count of ASVs",
        fill = "Primer Set"
    ) +
    theme_bw()

# Plot 4: Read abundance vs identity
p4 <- ggplot(all_blast, aes(x = Percent_Identity, y = Reads, color = Likely_Metazoan)) +
    geom_point(alpha = 0.6, size = 3) +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~PrimerSet) +
    labs(
        title = "ASV Abundance vs BLAST Identity",
        x = "Percent Identity",
        y = "Read Count (log scale)",
        color = "Likely Metazoan"
    ) +
    theme_bw()

# Save plots
ggsave(file.path(plots_dir, "hit_quality.pdf"), p1, width = 10, height = 6)
ggsave(file.path(plots_dir, "taxonomic_composition.pdf"), p2, width = 8, height = 6)
ggsave(file.path(plots_dir, "identity_distribution.pdf"), p3, width = 10, height = 8)
ggsave(file.path(plots_dir, "abundance_vs_identity.pdf"), p4, width = 10, height = 6)

# Combined plot
combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
        title = "BLAST Analysis of Unassigned ASVs",
        subtitle = sprintf("Analysis of %d unassigned ASVs", nrow(all_blast))
    )

ggsave(file.path(plots_dir, "blast_analysis_combined.pdf"),
       combined, width = 16, height = 12)

cat("\nPlots saved to:", plots_dir, "\n")

# Generate markdown report
report_file <- file.path(blast_dir, "BLAST_Analysis_Report.md")

sink(report_file)
cat("# BLAST Analysis Report: Unassigned ASVs\n\n")
cat(sprintf("**Date:** %s\n\n", Sys.Date()))
cat(sprintf("**Total ASVs Analyzed:** %d\n\n", nrow(all_blast)))

cat("## Summary Statistics\n\n")
cat("```\n")
print(summary_stats)
cat("```\n\n")

cat("## Key Findings\n\n")

# Calculate key metrics
pct_metazoan <- sum(all_blast$Likely_Metazoan) / nrow(all_blast) * 100
pct_coi <- sum(all_blast$Likely_COI) / nrow(all_blast) * 100
pct_good_hit <- sum(all_blast$Percent_Identity >= 90) / nrow(all_blast) * 100

cat(sprintf("- **%.1f%%** of unassigned ASVs are likely metazoan sequences\n", pct_metazoan))
cat(sprintf("- **%.1f%%** appear to be COI sequences based on hit descriptions\n", pct_coi))
cat(sprintf("- **%.1f%%** have ≥90%% identity to database sequences\n", pct_good_hit))
cat("\n")

cat("## Taxonomic Breakdown\n\n")
cat("```\n")
print(tax_breakdown)
cat("```\n\n")

cat("## Interpretation\n\n")

cat("### Possible reasons for unassigned status:\n\n")
cat("1. **Reference database gaps**: Sequences may be from species not well-represented in COInr\n")
cat("2. **Sequence quality issues**: Short or poor-quality sequences may lack sufficient overlap\n")
cat("3. **Novel diversity**: Potential new species or lineages\n")
cat("4. **Non-target amplification**: Contamination or non-COI amplicons\n")
cat("\n")

cat("### Recommendations:\n\n")

if (pct_metazoan > 50) {
    cat("- High proportion of metazoan hits suggests reference database could be expanded\n")
}

if (pct_coi > 70) {
    cat("- Most sequences appear to be COI - consider lowering confidence thresholds\n")
}

if (mean(all_blast$Percent_Identity) < 90) {
    cat("- Lower identity scores may indicate novel diversity or divergent lineages\n")
}

sink()

cat(sprintf("\nMarkdown report saved to: %s\n", report_file))

cat("\n================================================\n")
cat("Analysis Complete!\n")
cat("================================================\n\n")
