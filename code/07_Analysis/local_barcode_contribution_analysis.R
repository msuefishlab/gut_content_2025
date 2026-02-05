#!/usr/bin/env Rscript
# ==============================================================================
# Local Barcode Contribution Analysis
# ==============================================================================
# Purpose: Quantify the contribution of 429 locally-generated COI barcodes
#          (Schmidt lab) versus the existing COInr database to taxonomy
#          assignments in the metabarcoding pipeline.
#
# Author: Claude Code
# Date: 2026-01-30
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(qiime2R)
  library(scales)
})

# Set random seed for reproducibility
set.seed(42)

# ==============================================================================
# Section 1: Setup and Configuration
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("LOCAL BARCODE CONTRIBUTION ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Get root directory from environment
root <- Sys.getenv("root")
if (root == "") {
  stop("ERROR: 'root' environment variable not set. Source gut_contents.env first.")
}

# Define file paths
paths <- list(
  # Search results (QIIME2 artifacts)
  ps1_search = file.path(root, "output_data/05_Taxanomic_Assignment/primerset1_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.search_results.qza"),
  ps2_search = file.path(root, "output_data/05_Taxanomic_Assignment/primerset2_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.search_results.qza"),

  # Feature tables (for read abundance)
  ps1_table = file.path(root, "output_data/06_Generate_Output/primerset1_all_p985_table_filtd_NO_HOST.qza"),
  ps2_table = file.path(root, "output_data/06_Generate_Output/primerset2_all_p985_table_filtd_NO_HOST.qza"),

  # Taxonomy assignments
  ps1_taxonomy = file.path(root, "output_data/06_Generate_Output/primerset1_filtd_tax_dataframe_ALL.csv"),
  ps2_taxonomy = file.path(root, "output_data/06_Generate_Output/primerset2_filtd_tax_dataframe_ALL.csv"),

  # Local barcode reference
  local_barcodes = file.path(root, "output_data/04_Reference_Database/allschmidt429.tsv")
)

# Define output directory
output_dir <- file.path(root, "output_data/07_Analysis/local_barcode_contribution")
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Verify all input files exist
cat("Verifying input files...\n")
missing_files <- names(paths)[!sapply(paths, file.exists)]
if (length(missing_files) > 0) {
  stop(sprintf("ERROR: Missing input files:\n  %s",
               paste(missing_files, collapse = "\n  ")))
}
cat("✓ All input files found\n\n")

# ==============================================================================
# Section 2: Data Extraction and Classification Functions
# ==============================================================================

#' Extract VSEARCH results from QIIME2 artifact
#'
#' @param qza_path Path to search_results.qza file
#' @return Data frame with ASV-to-reference mappings
extract_search_results <- function(qza_path) {
  cat(sprintf("  Extracting: %s\n", basename(qza_path)))

  # Create temporary directory for extraction
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Unzip QZA file
  unzip(qza_path, exdir = temp_dir, overwrite = TRUE)

  # Find blast6.tsv file
  blast_file <- list.files(temp_dir,
                          pattern = "blast6\\.tsv$",
                          recursive = TRUE,
                          full.names = TRUE)[1]

  if (is.na(blast_file) || !file.exists(blast_file)) {
    stop(sprintf("ERROR: Could not find blast6.tsv in %s", qza_path))
  }

  # Read BLAST6 format results
  # Columns: query_id, reference_id, pident, length, mismatch, gapopen,
  #          qstart, qend, sstart, send, evalue, bitscore
  blast_data <- read_tsv(
    blast_file,
    col_names = c("ASV_ID", "Ref_ID", "Pident", "Length", "Mismatch",
                  "Gapopen", "Qstart", "Qend", "Sstart", "Send",
                  "Evalue", "Bitscore"),
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  ) %>%
    mutate(
      Pident = as.numeric(Pident),
      Length = as.integer(Length),
      Mismatch = as.integer(Mismatch),
      Evalue = as.numeric(Evalue),
      Bitscore = as.numeric(Bitscore)
    ) %>%
    # Remove no-match rows
    filter(Ref_ID != "*") %>%
    # Keep only top hit per ASV (already sorted by pident/bitscore in VSEARCH output)
    group_by(ASV_ID) %>%
    slice_head(n = 1) %>%
    ungroup()

  # Clean up
  unlink(temp_dir, recursive = TRUE)

  cat(sprintf("    Found %s ASVs with reference matches\n",
              comma(nrow(blast_data))))

  return(blast_data)
}

#' Classify reference source as Local or COInr
#'
#' @param ref_id Reference ID from BLAST results
#' @return Character vector: "Local" or "COInr"
classify_reference_source <- function(ref_id) {
  # Local barcodes end with "_all.fa" (e.g., "B001_1_all.fa", "3720_4_all.fa")
  # COInr references are numeric or GenBank accessions (e.g., "12633989", "MK642302_1")
  ifelse(str_ends(ref_id, "_all\\.fa"), "Local", "COInr")
}

# ==============================================================================
# Section 3: Data Loading
# ==============================================================================

cat("Loading datasets...\n")

# --- Load search results ---
cat("\n1. VSEARCH Results:\n")
ps1_search <- extract_search_results(paths$ps1_search)
ps2_search <- extract_search_results(paths$ps2_search)

# Classify sources
ps1_search <- ps1_search %>%
  mutate(Source = classify_reference_source(Ref_ID))

ps2_search <- ps2_search %>%
  mutate(Source = classify_reference_source(Ref_ID))

cat(sprintf("  PS1 - Local: %s (%.1f%%), COInr: %s (%.1f%%)\n",
            comma(sum(ps1_search$Source == "Local")),
            100 * mean(ps1_search$Source == "Local"),
            comma(sum(ps1_search$Source == "COInr")),
            100 * mean(ps1_search$Source == "COInr")))

cat(sprintf("  PS2 - Local: %s (%.1f%%), COInr: %s (%.1f%%)\n",
            comma(sum(ps2_search$Source == "Local")),
            100 * mean(ps2_search$Source == "Local"),
            comma(sum(ps2_search$Source == "COInr")),
            100 * mean(ps2_search$Source == "COInr")))

# --- Load feature tables ---
cat("\n2. Feature Tables (Read Abundance):\n")

ps1_table <- read_qza(paths$ps1_table)$data %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  mutate(Total_Reads = rowSums(select(., -ASV_ID))) %>%
  select(ASV_ID, Total_Reads)

ps2_table <- read_qza(paths$ps2_table)$data %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  mutate(Total_Reads = rowSums(select(., -ASV_ID))) %>%
  select(ASV_ID, Total_Reads)

cat(sprintf("  PS1: %s ASVs, %s total reads\n",
            comma(nrow(ps1_table)),
            comma(sum(ps1_table$Total_Reads))))

cat(sprintf("  PS2: %s ASVs, %s total reads\n",
            comma(nrow(ps2_table)),
            comma(sum(ps2_table$Total_Reads))))

# --- Load taxonomy ---
cat("\n3. Taxonomy Assignments:\n")

ps1_taxonomy <- read_csv(paths$ps1_taxonomy, show_col_types = FALSE) %>%
  rename(ASV_ID = `Feature ID`) %>%
  mutate(across(Kingdom:Species, ~str_remove(., "^[kpcofgs]__"))) %>%
  mutate(across(Kingdom:Species, ~str_remove(., "_[0-9]+$"))) %>%
  mutate(across(Kingdom:Species, ~na_if(., "")))

ps2_taxonomy <- read_csv(paths$ps2_taxonomy, show_col_types = FALSE) %>%
  rename(ASV_ID = `Feature ID`) %>%
  mutate(across(Kingdom:Species, ~str_remove(., "^[kpcofgs]__"))) %>%
  mutate(across(Kingdom:Species, ~str_remove(., "_[0-9]+$"))) %>%
  mutate(across(Kingdom:Species, ~na_if(., "")))

cat(sprintf("  PS1: %s taxonomic assignments\n", comma(nrow(ps1_taxonomy))))
cat(sprintf("  PS2: %s taxonomic assignments\n", comma(nrow(ps2_taxonomy))))

# --- Load local barcode reference ---
cat("\n4. Local Barcode Reference:\n")

local_barcodes <- read_tsv(paths$local_barcodes,
                          col_names = c("Barcode_ID", "Family_Ref", "Sequence"),
                          show_col_types = FALSE)

cat(sprintf("  %s local barcodes from Schmidt lab\n",
            comma(nrow(local_barcodes))))
cat(sprintf("  Families represented: %s\n",
            paste(sort(unique(local_barcodes$Family_Ref)), collapse = ", ")))

# ==============================================================================
# Section 4: Data Integration - Create Master Tables
# ==============================================================================

cat("\nCreating integrated master tables...\n")

create_master_table <- function(search_results, feature_table, taxonomy, primer_set) {

  master <- search_results %>%
    left_join(feature_table, by = "ASV_ID") %>%
    left_join(taxonomy, by = "ASV_ID") %>%
    mutate(
      Primer_Set = primer_set,
      Resolution = case_when(
        !is.na(Species) ~ "Species",
        !is.na(Genus) ~ "Genus",
        !is.na(Family) ~ "Family",
        !is.na(Order) ~ "Order",
        !is.na(Class) ~ "Class",
        !is.na(Phylum) ~ "Phylum",
        TRUE ~ "Kingdom"
      ),
      Resolution = factor(Resolution,
                         levels = c("Species", "Genus", "Family",
                                   "Order", "Class", "Phylum", "Kingdom"))
    )

  return(master)
}

ps1_master <- create_master_table(ps1_search, ps1_table, ps1_taxonomy, "PS1")
ps2_master <- create_master_table(ps2_search, ps2_table, ps2_taxonomy, "PS2")

# Combined table for cross-primer analyses
combined_master <- bind_rows(ps1_master, ps2_master)

cat(sprintf("  PS1 master table: %s rows\n", comma(nrow(ps1_master))))
cat(sprintf("  PS2 master table: %s rows\n", comma(nrow(ps2_master))))
cat(sprintf("  Combined table: %s rows\n", comma(nrow(combined_master))))

# Save ASV-level data (Table 1)
cat("\nSaving Table 01: ASV Database Hits...\n")
write_csv(
  combined_master %>%
    select(Primer_Set, ASV_ID, Ref_ID, Source, Pident, Total_Reads,
           Kingdom, Phylum, Class, Order, Family, Genus, Species, Resolution),
  file.path(output_dir, "tables/01_asv_database_hits.csv")
)

# ==============================================================================
# Section 5: Statistical Analyses
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("STATISTICAL ANALYSES\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# --- Analysis 1: Overall Contribution ---
cat("Analysis 1: Overall Contribution\n")
cat("─────────────────────────────────────────────────────────────────────\n")

overall_stats <- combined_master %>%
  group_by(Primer_Set, Source) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Primer_Set) %>%
  mutate(
    Percent_ASVs = 100 * N_ASVs / sum(N_ASVs),
    Percent_Reads = 100 * Total_Reads / sum(Total_Reads)
  ) %>%
  ungroup()

# Chi-square test: Are contributions significantly different?
ps1_contingency <- ps1_master %>%
  count(Source) %>%
  pull(n)

ps2_contingency <- ps2_master %>%
  count(Source) %>%
  pull(n)

chi_test_ps1 <- chisq.test(ps1_contingency)
chi_test_ps2 <- chisq.test(ps2_contingency)

cat("\nOverall Statistics:\n")
print(overall_stats %>%
        select(Primer_Set, Source, N_ASVs, Percent_ASVs,
               Total_Reads, Percent_Reads) %>%
        mutate(across(where(is.numeric) & !matches("Percent"), comma),
               across(matches("Percent"), ~sprintf("%.1f%%", .))))

cat(sprintf("\nChi-square test (PS1): χ² = %.2f, p = %.2e\n",
            chi_test_ps1$statistic, chi_test_ps1$p.value))
cat(sprintf("Chi-square test (PS2): χ² = %.2f, p = %.2e\n\n",
            chi_test_ps2$statistic, chi_test_ps2$p.value))

# Save Table 2
write_csv(overall_stats,
          file.path(output_dir, "tables/02_overall_contribution_summary.csv"))

# --- Analysis 2: Taxonomic Resolution Comparison ---
cat("Analysis 2: Taxonomic Resolution Comparison\n")
cat("─────────────────────────────────────────────────────────────────────\n")

resolution_stats <- combined_master %>%
  group_by(Primer_Set, Source, Resolution) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Primer_Set, Source) %>%
  mutate(
    Percent_ASVs = 100 * N_ASVs / sum(N_ASVs),
    Percent_Reads = 100 * Total_Reads / sum(Total_Reads)
  ) %>%
  ungroup()

# Chi-square: Does resolution depend on source?
resolution_contingency_ps1 <- ps1_master %>%
  count(Source, Resolution) %>%
  pivot_wider(names_from = Source, values_from = n, values_fill = 0) %>%
  select(-Resolution) %>%
  as.matrix()

resolution_contingency_ps2 <- ps2_master %>%
  count(Source, Resolution) %>%
  pivot_wider(names_from = Source, values_from = n, values_fill = 0) %>%
  select(-Resolution) %>%
  as.matrix()

chi_resolution_ps1 <- chisq.test(resolution_contingency_ps1)
chi_resolution_ps2 <- chisq.test(resolution_contingency_ps2)

cat("\nTaxonomic Resolution by Source:\n")
print(resolution_stats %>%
        arrange(Primer_Set, Source, Resolution) %>%
        select(Primer_Set, Source, Resolution, N_ASVs, Percent_Reads) %>%
        mutate(N_ASVs = comma(N_ASVs),
               Percent_Reads = sprintf("%.1f%%", Percent_Reads)))

cat(sprintf("\nChi-square test (PS1): χ² = %.2f, p = %.2e\n",
            chi_resolution_ps1$statistic, chi_resolution_ps1$p.value))
cat(sprintf("Chi-square test (PS2): χ² = %.2f, p = %.2e\n\n",
            chi_resolution_ps2$statistic, chi_resolution_ps2$p.value))

# Save Table 4
write_csv(resolution_stats,
          file.path(output_dir, "tables/04_taxonomic_resolution_comparison.csv"))

# --- Analysis 3: Family-Specific Breakdown ---
cat("Analysis 3: Family-Specific Breakdown\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Identify top families by read abundance
top_families <- combined_master %>%
  group_by(Family) %>%
  summarise(Total_Reads = sum(Total_Reads, na.rm = TRUE)) %>%
  arrange(desc(Total_Reads)) %>%
  filter(!is.na(Family)) %>%
  head(20) %>%
  pull(Family)

cat(sprintf("Analyzing top %d families by read abundance:\n", length(top_families)))
cat(sprintf("  %s\n\n", paste(top_families, collapse = ", ")))

family_stats <- combined_master %>%
  filter(Family %in% top_families) %>%
  group_by(Primer_Set, Family, Source) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Primer_Set, Family) %>%
  mutate(
    Percent_ASVs_in_Family = 100 * N_ASVs / sum(N_ASVs),
    Percent_Reads_in_Family = 100 * Total_Reads / sum(Total_Reads),
    Total_Family_Reads = sum(Total_Reads)
  ) %>%
  ungroup() %>%
  arrange(Primer_Set, desc(Total_Family_Reads), Family, Source)

cat("Family-Level Contribution Summary (top 10):\n")
print(family_stats %>%
        head(20) %>%
        select(Primer_Set, Family, Source, N_ASVs,
               Percent_Reads_in_Family, Total_Family_Reads) %>%
        mutate(N_ASVs = comma(N_ASVs),
               Percent_Reads_in_Family = sprintf("%.1f%%", Percent_Reads_in_Family),
               Total_Family_Reads = comma(Total_Family_Reads)))

cat("\n")

# Save Table 3
write_csv(family_stats,
          file.path(output_dir, "tables/03_family_contribution_summary.csv"))

# --- Analysis 4: Database Gap Analysis ---
cat("Analysis 4: Database Gap Analysis\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Families where local barcodes are the ONLY source
local_only <- combined_master %>%
  filter(!is.na(Family)) %>%
  group_by(Primer_Set, Family) %>%
  summarise(
    Has_Local = any(Source == "Local"),
    Has_COInr = any(Source == "COInr"),
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Has_Local & !Has_COInr) %>%
  arrange(Primer_Set, desc(Total_Reads))

cat("Families with ONLY local barcode coverage:\n")
if (nrow(local_only) > 0) {
  print(local_only %>%
          select(Primer_Set, Family, N_ASVs, Total_Reads) %>%
          mutate(N_ASVs = comma(N_ASVs),
                 Total_Reads = comma(Total_Reads)))
} else {
  cat("  None found\n")
}
cat("\n")

# Families where local barcodes are dominant (>50% of reads)
local_dominant <- combined_master %>%
  filter(!is.na(Family)) %>%
  group_by(Primer_Set, Family, Source) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Primer_Set, Family) %>%
  mutate(
    Percent_Reads = 100 * Total_Reads / sum(Total_Reads),
    Total_Family_Reads = sum(Total_Reads)
  ) %>%
  ungroup() %>%
  filter(Source == "Local", Percent_Reads > 50) %>%
  arrange(Primer_Set, desc(Total_Family_Reads))

cat("Families with >50% local barcode contribution:\n")
if (nrow(local_dominant) > 0) {
  print(local_dominant %>%
          select(Primer_Set, Family, N_ASVs, Percent_Reads, Total_Family_Reads) %>%
          mutate(N_ASVs = comma(N_ASVs),
                 Percent_Reads = sprintf("%.1f%%", Percent_Reads),
                 Total_Family_Reads = comma(Total_Family_Reads)))
} else {
  cat("  None found\n")
}
cat("\n")

# Save Table 5
write_csv(bind_rows(
  local_only %>% mutate(Category = "Local Only"),
  local_dominant %>% mutate(Category = "Local Dominant (>50%)")
),
file.path(output_dir, "tables/05_local_only_families.csv"))

# --- Analysis 5: Cross-Primer Comparison ---
cat("Analysis 5: Cross-Primer Comparison (PS1 vs PS2)\n")
cat("─────────────────────────────────────────────────────────────────────\n")

cross_primer <- combined_master %>%
  group_by(Primer_Set, Source) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Primer_Set) %>%
  mutate(
    Percent_ASVs = 100 * N_ASVs / sum(N_ASVs),
    Percent_Reads = 100 * Total_Reads / sum(Total_Reads)
  ) %>%
  ungroup()

cat("Cross-Primer Consistency:\n")
print(cross_primer %>%
        select(Primer_Set, Source, Percent_ASVs, Percent_Reads) %>%
        mutate(Percent_ASVs = sprintf("%.1f%%", Percent_ASVs),
               Percent_Reads = sprintf("%.1f%%", Percent_Reads)))

cat("\n")

# Save Table 6
write_csv(cross_primer,
          file.path(output_dir, "tables/06_ps1_vs_ps2_comparison.csv"))

# Save statistical tests summary
statistical_tests <- tibble(
  Test = c("Overall Contribution (PS1)", "Overall Contribution (PS2)",
           "Taxonomic Resolution (PS1)", "Taxonomic Resolution (PS2)"),
  Chi_Square = c(chi_test_ps1$statistic, chi_test_ps2$statistic,
                 chi_resolution_ps1$statistic, chi_resolution_ps2$statistic),
  P_Value = c(chi_test_ps1$p.value, chi_test_ps2$p.value,
              chi_resolution_ps1$p.value, chi_resolution_ps2$p.value),
  Significant = ifelse(P_Value < 0.05, "Yes", "No")
)

write_csv(statistical_tests,
          file.path(output_dir, "tables/statistical_tests_summary.csv"))

# ==============================================================================
# Section 6: Visualizations
# ==============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("CREATING VISUALIZATIONS\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Define colorblind-friendly palette
colors_source <- c("Local" = "#E69F00", "COInr" = "#56B4E9")

# --- Plot 1: Overall Contribution (Stacked Bar) ---
cat("Plot 1: Overall Contribution (ASVs and Reads)...\n")

plot_data_1 <- overall_stats %>%
  select(Primer_Set, Source, Percent_ASVs, Percent_Reads) %>%
  pivot_longer(cols = c(Percent_ASVs, Percent_Reads),
               names_to = "Metric",
               values_to = "Percentage") %>%
  mutate(Metric = recode(Metric,
                        "Percent_ASVs" = "% ASVs",
                        "Percent_Reads" = "% Reads"))

p1 <- ggplot(plot_data_1, aes(x = Metric, y = Percentage, fill = Source)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  facet_wrap(~Primer_Set, ncol = 2) +
  scale_fill_manual(values = colors_source) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "Overall Contribution of Local Barcodes vs COInr Database",
       subtitle = "Proportion of ASVs and reads assigned to each database source",
       x = NULL,
       y = "Percentage (%)",
       fill = "Database Source") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "figures/01_overall_contribution.pdf"),
       p1, width = 8, height = 6)

# --- Plot 2: Family-Specific Breakdown (Dodged Bar) ---
cat("Plot 2: Family-Specific Breakdown...\n")

plot_data_2 <- family_stats %>%
  filter(Family %in% top_families[1:10]) %>%
  mutate(Family = fct_reorder(Family, Total_Family_Reads, .desc = TRUE))

p2 <- ggplot(plot_data_2, aes(x = Family, y = Total_Reads, fill = Source)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~Primer_Set, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = colors_source) +
  scale_y_log10(labels = comma) +
  labs(title = "Family-Specific Database Contribution",
       subtitle = "Total reads by database source for top 10 families",
       x = "Family",
       y = "Total Reads (log scale)",
       fill = "Database Source") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "figures/02_family_breakdown.pdf"),
       p2, width = 10, height = 8)

# --- Plot 3: Taxonomic Resolution Comparison (Dodged Bar) ---
cat("Plot 3: Taxonomic Resolution Comparison...\n")

plot_data_3 <- resolution_stats %>%
  group_by(Primer_Set, Resolution) %>%
  mutate(Total_Resolution_Reads = sum(Total_Reads)) %>%
  ungroup()

p3 <- ggplot(plot_data_3, aes(x = Resolution, y = Percent_Reads, fill = Source)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Percent_Reads)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3) +
  facet_wrap(~Primer_Set, ncol = 2) +
  scale_fill_manual(values = colors_source) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Taxonomic Resolution by Database Source",
       subtitle = "Percentage of reads assigned to each taxonomic level",
       x = "Taxonomic Resolution",
       y = "% of Reads",
       fill = "Database Source") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "figures/03_resolution_comparison.pdf"),
       p3, width = 10, height = 6)

# --- Plot 4: Cross-Primer Consistency (Comparison) ---
cat("Plot 4: Cross-Primer Consistency (PS1 vs PS2)...\n")

plot_data_4 <- cross_primer %>%
  select(Primer_Set, Source, Percent_Reads) %>%
  pivot_wider(names_from = Primer_Set, values_from = Percent_Reads)

p4 <- ggplot(cross_primer, aes(x = Primer_Set, y = Percent_Reads, fill = Source)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Reads)),
            position = position_dodge(width = 0.6),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = colors_source) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Cross-Primer Set Consistency",
       subtitle = "Local barcode contribution across primer sets (by reads)",
       x = "Primer Set",
       y = "% of Reads",
       fill = "Database Source") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "figures/04_ps1_vs_ps2_patterns.pdf"),
       p4, width = 7, height = 6)

cat("✓ All figures saved\n\n")

# ==============================================================================
# Section 7: Generate Markdown Report
# ==============================================================================

cat("Generating comprehensive report...\n")

# Calculate key statistics for report
ps1_local_reads_pct <- overall_stats %>%
  filter(Primer_Set == "PS1", Source == "Local") %>%
  pull(Percent_Reads)

ps2_local_reads_pct <- overall_stats %>%
  filter(Primer_Set == "PS2", Source == "Local") %>%
  pull(Percent_Reads)

total_asvs <- nrow(combined_master)
local_asvs_pct <- 100 * sum(combined_master$Source == "Local") / total_asvs

# Create report
report_lines <- c(
  "# Local Barcode Contribution Analysis",
  "",
  sprintf("**Analysis Date:** %s", Sys.Date()),
  sprintf("**Total ASVs Analyzed:** %s", comma(total_asvs)),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  "This analysis quantifies the contribution of 429 locally-generated COI barcodes from the Schmidt lab versus the existing COInr public database to taxonomy assignments in our gut content metabarcoding pipeline.",
  "",
  "### Key Findings",
  "",
  sprintf("- **Overall ASV Contribution:** %.1f%% of ASVs matched local barcodes", local_asvs_pct),
  sprintf("- **Read Abundance (PS1):** %.1f%% of reads were assigned using local barcodes", ps1_local_reads_pct),
  sprintf("- **Read Abundance (PS2):** %.1f%% of reads were assigned using local barcodes", ps2_local_reads_pct),
  sprintf("- **Statistical Significance:** Local vs COInr contributions are %s (p < 0.05)",
          ifelse(chi_test_ps1$p.value < 0.05, "significantly different", "not significantly different")),
  "",
  if (nrow(local_only) > 0) {
    sprintf("- **Database Gaps Filled:** %d families had coverage ONLY from local barcodes",
            length(unique(local_only$Family)))
  } else {
    "- **Database Gaps Filled:** No families relied exclusively on local barcodes"
  },
  "",
  "### Conclusion",
  "",
  if (ps1_local_reads_pct > 20 || ps2_local_reads_pct > 20) {
    sprintf("The investment in local barcode generation provided **substantial value**, with local barcodes contributing %.1f%% (PS1) and %.1f%% (PS2) of total read assignments. This demonstrates that locally-generated barcodes filled critical gaps in the public database for regionally-relevant prey taxa.",
            ps1_local_reads_pct, ps2_local_reads_pct)
  } else if (ps1_local_reads_pct > 5 || ps2_local_reads_pct > 5) {
    sprintf("The local barcode library provided **moderate value**, contributing %.1f%% (PS1) and %.1f%% (PS2) of read assignments. While not the dominant source, local barcodes improved taxonomic coverage for specific families.",
            ps1_local_reads_pct, ps2_local_reads_pct)
  } else {
    sprintf("The local barcode contribution was **limited** (%.1f%% PS1, %.1f%% PS2 of reads), suggesting the COInr database already had good coverage for most regionally-relevant taxa. However, local barcodes may have filled specific gaps for rare or underrepresented families.",
            ps1_local_reads_pct, ps2_local_reads_pct)
  },
  "",
  "---",
  "",
  "## 1. Overall Contribution",
  "",
  "### Summary Statistics",
  "",
  "| Primer Set | Database | ASVs | % ASVs | Reads | % Reads |",
  "|------------|----------|------|--------|-------|---------|",
  paste0(apply(overall_stats %>%
                 select(Primer_Set, Source, N_ASVs, Percent_ASVs, Total_Reads, Percent_Reads) %>%
                 mutate(N_ASVs = comma(N_ASVs),
                        Percent_ASVs = sprintf("%.1f%%", Percent_ASVs),
                        Total_Reads = comma(Total_Reads),
                        Percent_Reads = sprintf("%.1f%%", Percent_Reads)),
               1, function(x) sprintf("| %s | %s | %s | %s | %s | %s |", x[1], x[2], x[3], x[4], x[5], x[6]))),
  "",
  "### Statistical Test",
  "",
  sprintf("- **Chi-square test (PS1):** χ² = %.2f, p = %.2e",
          chi_test_ps1$statistic, chi_test_ps1$p.value),
  sprintf("- **Chi-square test (PS2):** χ² = %.2f, p = %.2e",
          chi_test_ps2$statistic, chi_test_ps2$p.value),
  sprintf("- **Interpretation:** The difference between local and COInr contributions is %s.",
          ifelse(chi_test_ps1$p.value < 0.05 | chi_test_ps2$p.value < 0.05,
                 "**statistically significant**",
                 "not statistically significant")),
  "",
  "### Visualization",
  "",
  "![Overall Contribution](figures/01_overall_contribution.pdf)",
  "",
  "---",
  "",
  "## 2. Taxonomic Resolution",
  "",
  "### Resolution by Database Source (% Reads)",
  "",
  "| Primer Set | Source | Species | Genus | Family | Order | Class | Phylum |",
  "|------------|--------|---------|-------|--------|-------|-------|--------|",
  paste0(apply(resolution_stats %>%
                 select(Primer_Set, Source, Resolution, Percent_Reads) %>%
                 pivot_wider(names_from = Resolution, values_from = Percent_Reads, values_fill = 0) %>%
                 mutate(across(where(is.numeric), ~sprintf("%.1f%%", .))),
               1, function(x) {
                 sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                        x["Primer_Set"], x["Source"],
                        ifelse("Species" %in% names(x), x["Species"], "0.0%"),
                        ifelse("Genus" %in% names(x), x["Genus"], "0.0%"),
                        ifelse("Family" %in% names(x), x["Family"], "0.0%"),
                        ifelse("Order" %in% names(x), x["Order"], "0.0%"),
                        ifelse("Class" %in% names(x), x["Class"], "0.0%"),
                        ifelse("Phylum" %in% names(x), x["Phylum"], "0.0%"))
               })),
  "",
  "### Statistical Test",
  "",
  sprintf("- **Chi-square test (PS1):** χ² = %.2f, p = %.2e",
          chi_resolution_ps1$statistic, chi_resolution_ps1$p.value),
  sprintf("- **Chi-square test (PS2):** χ² = %.2f, p = %.2e",
          chi_resolution_ps2$statistic, chi_resolution_ps2$p.value),
  "",
  "### Visualization",
  "",
  "![Taxonomic Resolution](figures/03_resolution_comparison.pdf)",
  "",
  "---",
  "",
  "## 3. Family-Specific Analysis",
  "",
  sprintf("Analysis focused on the top %d families by read abundance:", length(top_families)),
  "",
  paste0("- ", top_families),
  "",
  "### Top Families by Total Reads",
  "",
  "| Primer Set | Family | Source | ASVs | % Reads in Family | Total Family Reads |",
  "|------------|--------|--------|------|-------------------|--------------------|",
  paste0(apply(family_stats %>%
                 head(20) %>%
                 select(Primer_Set, Family, Source, N_ASVs,
                       Percent_Reads_in_Family, Total_Family_Reads) %>%
                 mutate(N_ASVs = comma(N_ASVs),
                        Percent_Reads_in_Family = sprintf("%.1f%%", Percent_Reads_in_Family),
                        Total_Family_Reads = comma(Total_Family_Reads)),
               1, function(x) sprintf("| %s | %s | %s | %s | %s | %s |",
                                     x[1], x[2], x[3], x[4], x[5], x[6]))),
  "",
  "### Visualization",
  "",
  "![Family Breakdown](figures/02_family_breakdown.pdf)",
  "",
  "---",
  "",
  "## 4. Database Gaps Filled by Local Barcodes",
  "",
  "### Families with ONLY Local Barcode Coverage",
  "",
  if (nrow(local_only) > 0) {
    c(
      "| Primer Set | Family | ASVs | Total Reads |",
      "|------------|--------|------|-------------|",
      paste0(apply(local_only %>%
                     select(Primer_Set, Family, N_ASVs, Total_Reads) %>%
                     mutate(N_ASVs = comma(N_ASVs),
                            Total_Reads = comma(Total_Reads)),
                   1, function(x) sprintf("| %s | %s | %s | %s |", x[1], x[2], x[3], x[4])))
    )
  } else {
    "No families were assigned exclusively using local barcodes. The COInr database had at least some representation for all detected families."
  },
  "",
  "### Families with >50% Local Contribution",
  "",
  if (nrow(local_dominant) > 0) {
    c(
      "| Primer Set | Family | ASVs | % Reads | Total Family Reads |",
      "|------------|--------|------|---------|-------------------|",
      paste0(apply(local_dominant %>%
                     select(Primer_Set, Family, N_ASVs, Percent_Reads, Total_Family_Reads) %>%
                     mutate(N_ASVs = comma(N_ASVs),
                            Percent_Reads = sprintf("%.1f%%", Percent_Reads),
                            Total_Family_Reads = comma(Total_Family_Reads)),
                   1, function(x) sprintf("| %s | %s | %s | %s | %s |", x[1], x[2], x[3], x[4], x[5])))
    )
  } else {
    "No families had >50% of reads assigned via local barcodes."
  },
  "",
  "---",
  "",
  "## 5. Cross-Primer Consistency",
  "",
  "### PS1 vs PS2 Comparison",
  "",
  "| Primer Set | Source | % ASVs | % Reads |",
  "|------------|--------|--------|---------|",
  paste0(apply(cross_primer %>%
                 select(Primer_Set, Source, Percent_ASVs, Percent_Reads) %>%
                 mutate(Percent_ASVs = sprintf("%.1f%%", Percent_ASVs),
                        Percent_Reads = sprintf("%.1f%%", Percent_Reads)),
               1, function(x) sprintf("| %s | %s | %s | %s |", x[1], x[2], x[3], x[4]))),
  "",
  "### Interpretation",
  "",
  sprintf("The local barcode contribution was %s across primer sets, with PS1 showing %.1f%% and PS2 showing %.1f%% of reads assigned via local barcodes.",
          ifelse(abs(ps1_local_reads_pct - ps2_local_reads_pct) < 5, "**highly consistent**", "**moderately variable**"),
          ps1_local_reads_pct, ps2_local_reads_pct),
  "",
  "### Visualization",
  "",
  "![Cross-Primer Patterns](figures/04_ps1_vs_ps2_patterns.pdf)",
  "",
  "---",
  "",
  "## Methods",
  "",
  "### Pipeline Parameters",
  "",
  "- **Search algorithm:** VSEARCH",
  "- **Identity threshold:** 95%",
  "- **Coverage threshold:** 94%",
  "- **Reference database:** COInr Metazoa + 429 Schmidt lab barcodes (Leray-trimmed)",
  "",
  "### Statistical Tests",
  "",
  "- **Chi-square tests** for overall contribution and taxonomic resolution",
  "- **Significance threshold:** p < 0.05",
  "",
  "### Data Processing",
  "",
  "1. Extracted top VSEARCH hit per ASV from `search_results.qza`",
  "2. Classified references as 'Local' (ending in `_all.fa`) or 'COInr' (all others)",
  "3. Integrated with feature tables for read abundance",
  "4. Integrated with taxonomy assignments for resolution analysis",
  "",
  "### Software",
  "",
  "- QIIME2 (for metabarcoding pipeline)",
  "- R 4.x (tidyverse, qiime2R, scales)",
  "",
  "---",
  "",
  "## Output Files",
  "",
  "### Tables",
  "",
  "1. `01_asv_database_hits.csv` - Complete ASV-to-reference mappings with source classification",
  "2. `02_overall_contribution_summary.csv` - Overall statistics by primer set and source",
  "3. `03_family_contribution_summary.csv` - Family-level breakdown of contributions",
  "4. `04_taxonomic_resolution_comparison.csv` - Resolution levels by source",
  "5. `05_local_only_families.csv` - Families with exclusive or dominant local coverage",
  "6. `06_ps1_vs_ps2_comparison.csv` - Cross-primer set comparison",
  "7. `statistical_tests_summary.csv` - Chi-square test results",
  "",
  "### Figures",
  "",
  "1. `01_overall_contribution.pdf` - Stacked bar chart of ASV/read percentages",
  "2. `02_family_breakdown.pdf` - Family-specific source contributions",
  "3. `03_resolution_comparison.pdf` - Taxonomic resolution by source",
  "4. `04_ps1_vs_ps2_patterns.pdf` - Cross-primer consistency",
  "",
  "---",
  "",
  sprintf("**Report generated:** %s", Sys.time()),
  sprintf("**Analysis script:** `code/07_Analysis/local_barcode_contribution_analysis.R`"),
  ""
)

# Write report
writeLines(report_lines,
          file.path(output_dir, "local_barcode_contribution_report.md"))

cat("✓ Report saved\n\n")

# ==============================================================================
# Final Summary
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat("Output Directory:\n")
cat(sprintf("  %s\n\n", output_dir))

cat("Generated Files:\n")
cat("  Tables (7):\n")
cat("    ✓ 01_asv_database_hits.csv\n")
cat("    ✓ 02_overall_contribution_summary.csv\n")
cat("    ✓ 03_family_contribution_summary.csv\n")
cat("    ✓ 04_taxonomic_resolution_comparison.csv\n")
cat("    ✓ 05_local_only_families.csv\n")
cat("    ✓ 06_ps1_vs_ps2_comparison.csv\n")
cat("    ✓ statistical_tests_summary.csv\n\n")

cat("  Figures (4):\n")
cat("    ✓ 01_overall_contribution.pdf\n")
cat("    ✓ 02_family_breakdown.pdf\n")
cat("    ✓ 03_resolution_comparison.pdf\n")
cat("    ✓ 04_ps1_vs_ps2_patterns.pdf\n\n")

cat("  Report:\n")
cat("    ✓ local_barcode_contribution_report.md\n\n")

cat("Review the comprehensive report at:\n")
cat(sprintf("  %s\n\n", file.path(output_dir, "local_barcode_contribution_report.md")))

cat("═══════════════════════════════════════════════════════════════════════\n")
