#!/usr/bin/env Rscript

# ==============================================================================
# Taxonomy Assignment Rate Analysis: PS1 vs PS2
# ==============================================================================
# Purpose: Investigate whether differential taxonomy assignment rates between
#          COI primer sets could explain observed differences in prey composition
#
# Author: Claude Code
# Date: 2026-01-30
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(qiime2R)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Setup and Configuration
# ==============================================================================

# Get repository root from environment variable
root <- Sys.getenv("root")
if (root == "") {
  stop("Error: 'root' environment variable not set. Source gut_contents.env first.")
}

# Define paths
output_base <- file.path(root, "output_data/07_Analysis/assignment_diagnostics")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "sequences"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "figures"), recursive = TRUE, showWarnings = FALSE)

# Input file paths
ps1_tax_file <- file.path(root, "output_data/06_Generate_Output/primerset1_filtd_tax_dataframe_ALL.csv")
ps2_tax_file <- file.path(root, "output_data/06_Generate_Output/primerset2_filtd_tax_dataframe_ALL.csv")
ps1_table_file <- file.path(root, "output_data/06_Generate_Output/primerset1_all_p985_table_filtd_NO_HOST.qza")
ps2_table_file <- file.path(root, "output_data/06_Generate_Output/primerset2_all_p985_table_filtd_NO_HOST.qza")
ps1_seqs_file <- file.path(root, "output_data/06_Generate_Output/primerset1_all_p985_seqs_Filtd.qza")
ps2_seqs_file <- file.path(root, "output_data/06_Generate_Output/primerset2_all_p985_seqs_Filtd.qza")
ps1_metadata_file <- file.path(root, "input_data/metadata/primerset1_metadata_final.tsv")
ps2_metadata_file <- file.path(root, "input_data/metadata/primerset2_metadata_final.tsv")

cat("Starting taxonomy assignment diagnostics analysis...\n")
cat("Output directory:", output_base, "\n\n")

# ==============================================================================
# Data Loading Functions
# ==============================================================================

load_taxonomy <- function(tax_file) {
  cat("Loading taxonomy from:", basename(tax_file), "\n")
  tax <- read_csv(tax_file, show_col_types = FALSE)

  # Rename column with space to use period
  tax <- tax %>%
    rename(Feature.ID = `Feature ID`)

  # Clean taxonomy columns - remove prefixes and taxids, trim whitespace
  tax_parsed <- tax %>%
    mutate(across(Kingdom:Species, ~str_trim(.))) %>%
    mutate(across(Kingdom:Species, ~str_remove(., "^[kpcofgs]__"))) %>%
    mutate(across(Kingdom:Species, ~str_remove(., "_[0-9]+$"))) %>%
    mutate(across(Kingdom:Species, ~na_if(., ""))) %>%
    mutate(across(Kingdom:Species, ~na_if(., "NA")))

  cat("  Loaded", nrow(tax_parsed), "ASVs\n")
  return(tax_parsed)
}

load_feature_table <- function(qza_file) {
  cat("Loading feature table from:", basename(qza_file), "\n")
  ft <- read_qza(qza_file)$data
  ft_df <- as.data.frame(ft) %>%
    rownames_to_column("ASV") %>%
    as_tibble()
  cat("  Loaded", nrow(ft_df), "ASVs x", ncol(ft_df)-1, "samples\n")
  return(ft_df)
}

load_sequences <- function(qza_file) {
  cat("Loading sequences from:", basename(qza_file), "\n")
  seqs <- read_qza(qza_file)$data
  cat("  Loaded", length(seqs), "sequences\n")
  return(seqs)
}

# ==============================================================================
# Load All Data
# ==============================================================================

cat("\n=== Loading Data ===\n")
ps1_tax <- load_taxonomy(ps1_tax_file)
ps2_tax <- load_taxonomy(ps2_tax_file)
ps1_table <- load_feature_table(ps1_table_file)
ps2_table <- load_feature_table(ps2_table_file)
ps1_seqs <- load_sequences(ps1_seqs_file)
ps2_seqs <- load_sequences(ps2_seqs_file)
ps1_metadata <- read_tsv(ps1_metadata_file, show_col_types = FALSE)
ps2_metadata <- read_tsv(ps2_metadata_file, show_col_types = FALSE)

# Validation checks
cat("\n=== Validation Checks ===\n")
cat("PS1 ASVs: Taxonomy =", nrow(ps1_tax), ", Table =", nrow(ps1_table), ", Seqs =", length(ps1_seqs), "\n")
cat("PS2 ASVs: Taxonomy =", nrow(ps2_tax), ", Table =", nrow(ps2_table), ", Seqs =", length(ps2_seqs), "\n")

# Initialize results storage
stats_tests <- list()

# ==============================================================================
# Analysis 1: ASV-Level Assignment Rates
# ==============================================================================

cat("\n=== Analysis 1: ASV-Level Assignment Rates ===\n")

calc_assignment_rates <- function(tax_df, primer_set) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  rates <- tibble(
    Rank = ranks,
    PrimerSet = primer_set,
    Total_ASVs = nrow(tax_df),
    Assigned_ASVs = sapply(ranks, function(r) sum(!is.na(tax_df[[r]]))),
    Percent_Assigned = round(100 * Assigned_ASVs / Total_ASVs, 2)
  )

  return(rates)
}

ps1_rates <- calc_assignment_rates(ps1_tax, "PS1")
ps2_rates <- calc_assignment_rates(ps2_tax, "PS2")
asv_rates <- bind_rows(ps1_rates, ps2_rates)

write_csv(asv_rates, file.path(output_base, "tables/01_asv_assignment_rates.csv"))
cat("Saved: tables/01_asv_assignment_rates.csv\n")

# Statistical test: Chi-square at Family level
ps1_family <- sum(!is.na(ps1_tax$Family))
ps1_no_family <- sum(is.na(ps1_tax$Family))
ps2_family <- sum(!is.na(ps2_tax$Family))
ps2_no_family <- sum(is.na(ps2_tax$Family))

chi_test <- chisq.test(matrix(c(ps1_family, ps1_no_family, ps2_family, ps2_no_family), nrow=2))
stats_tests$asv_assignment <- tibble(
  Test = "Chi-square: ASV Family Assignment (PS1 vs PS2)",
  Statistic = chi_test$statistic,
  P_value = chi_test$p.value,
  Effect_Size = sqrt(chi_test$statistic / sum(c(ps1_family, ps1_no_family, ps2_family, ps2_no_family)))
)

# Plot
p1 <- ggplot(asv_rates, aes(x = Rank, y = Percent_Assigned, fill = PrimerSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Assigned)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "ASV-Level Assignment Rates by Taxonomic Rank",
       subtitle = sprintf("PS1: n=%d ASVs | PS2: n=%d ASVs", nrow(ps1_tax), nrow(ps2_tax)),
       x = "Taxonomic Rank", y = "% ASVs Assigned", fill = "Primer Set") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_base, "figures/01_assignment_rates_by_rank.pdf"), p1, width = 8, height = 6)
cat("Saved: figures/01_assignment_rates_by_rank.pdf\n")

# ==============================================================================
# Analysis 2: Read-Weighted Assignment Rates
# ==============================================================================

cat("\n=== Analysis 2: Read-Weighted Assignment Rates ===\n")

calc_read_weighted_rates <- function(tax_df, table_df) {
  # Merge taxonomy with abundance
  merged <- tax_df %>%
    select(Feature.ID, Family) %>%
    inner_join(table_df, by = c("Feature.ID" = "ASV"))

  # Global rate
  total_reads <- sum(select(merged, -Feature.ID, -Family))
  family_reads <- merged %>%
    filter(!is.na(Family)) %>%
    select(-Feature.ID, -Family) %>%
    sum()

  global_rate <- 100 * family_reads / total_reads

  # Per-sample rates
  sample_rates <- merged %>%
    select(-Feature.ID) %>%
    pivot_longer(-Family, names_to = "Sample", values_to = "Reads") %>%
    group_by(Sample) %>%
    summarise(
      Total_Reads = sum(Reads),
      Family_Reads = sum(Reads[!is.na(Family)]),
      Percent_Assigned = 100 * Family_Reads / Total_Reads,
      .groups = "drop"
    )

  return(list(global = global_rate, per_sample = sample_rates))
}

ps1_read_rates <- calc_read_weighted_rates(ps1_tax, ps1_table)
ps2_read_rates <- calc_read_weighted_rates(ps2_tax, ps2_table)

# Global rates
global_rates <- tibble(
  PrimerSet = c("PS1", "PS2"),
  Percent_Reads_Assigned_Family = c(ps1_read_rates$global, ps2_read_rates$global)
)
write_csv(global_rates, file.path(output_base, "tables/02_read_weighted_assignment_global.csv"))
cat("Saved: tables/02_read_weighted_assignment_global.csv\n")

# Per-sample rates
per_sample_rates <- bind_rows(
  ps1_read_rates$per_sample %>% mutate(PrimerSet = "PS1"),
  ps2_read_rates$per_sample %>% mutate(PrimerSet = "PS2")
)
write_csv(per_sample_rates, file.path(output_base, "tables/02_read_weighted_assignment_per_sample.csv"))
cat("Saved: tables/02_read_weighted_assignment_per_sample.csv\n")

# Statistical test: Wilcoxon rank-sum
wilcox_test <- wilcox.test(
  ps1_read_rates$per_sample$Percent_Assigned,
  ps2_read_rates$per_sample$Percent_Assigned
)

stats_tests$read_weighted <- tibble(
  Test = "Wilcoxon: Per-Sample Family Assignment (PS1 vs PS2)",
  Statistic = wilcox_test$statistic,
  P_value = wilcox_test$p.value,
  Effect_Size = NA_real_
)

# Plot
p2 <- ggplot(per_sample_rates, aes(x = PrimerSet, y = Percent_Assigned, fill = PrimerSet)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "Per-Sample Read-Weighted Family Assignment Rates",
       subtitle = sprintf("PS1: %.1f%% | PS2: %.1f%% | p = %.3f",
                         ps1_read_rates$global, ps2_read_rates$global, wilcox_test$p.value),
       x = "Primer Set", y = "% Reads Assigned to Family") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(output_base, "figures/02_read_weighted_boxplot.pdf"), p2, width = 6, height = 6)
cat("Saved: figures/02_read_weighted_boxplot.pdf\n")

# ==============================================================================
# Analysis 3: Taxon-Specific Analysis
# ==============================================================================

cat("\n=== Analysis 3: Taxon-Specific Analysis ===\n")

calc_family_stats <- function(tax_df, table_df, primer_set) {
  # Merge taxonomy with abundance
  merged <- tax_df %>%
    select(Feature.ID, Kingdom, Phylum, Class, Order, Family) %>%
    inner_join(table_df, by = c("Feature.ID" = "ASV"))

  # Calculate total abundance per ASV
  merged <- merged %>%
    mutate(Total_Reads = rowSums(select(., -Feature.ID, -Kingdom, -Phylum, -Class, -Order, -Family)))

  # Family-level stats
  family_stats <- merged %>%
    filter(!is.na(Family)) %>%
    group_by(Family) %>%
    summarise(
      N_ASVs = n(),
      Total_Reads = sum(Total_Reads),
      .groups = "drop"
    ) %>%
    mutate(
      PrimerSet = primer_set,
      Rel_Abundance = 100 * Total_Reads / sum(Total_Reads)
    ) %>%
    arrange(desc(Total_Reads))

  # Order-only assignments (not assigned to Family)
  order_only <- merged %>%
    filter(!is.na(Order) & is.na(Family)) %>%
    group_by(Order) %>%
    summarise(
      N_ASVs = n(),
      Total_Reads = sum(Total_Reads),
      .groups = "drop"
    ) %>%
    mutate(
      PrimerSet = primer_set,
      Rel_Abundance = 100 * Total_Reads / sum(merged$Total_Reads)
    ) %>%
    arrange(desc(Total_Reads))

  return(list(family = family_stats, order_only = order_only))
}

ps1_family_stats <- calc_family_stats(ps1_tax, ps1_table, "PS1")
ps2_family_stats <- calc_family_stats(ps2_tax, ps2_table, "PS2")

# Combine family stats
family_stats_combined <- bind_rows(ps1_family_stats$family, ps2_family_stats$family)
write_csv(family_stats_combined, file.path(output_base, "tables/03_family_specific_stats.csv"))
cat("Saved: tables/03_family_specific_stats.csv\n")

# Combine order-only stats
order_only_combined <- bind_rows(ps1_family_stats$order_only, ps2_family_stats$order_only)
write_csv(order_only_combined, file.path(output_base, "tables/03_order_only_assignments.csv"))
cat("Saved: tables/03_order_only_assignments.csv\n")

# Plot: Top families comparison
top_families <- family_stats_combined %>%
  group_by(Family) %>%
  summarise(Max_Reads = max(Total_Reads), .groups = "drop") %>%
  slice_max(Max_Reads, n = 10) %>%
  pull(Family)

p3 <- family_stats_combined %>%
  filter(Family %in% top_families) %>%
  ggplot(aes(x = reorder(Family, -Total_Reads), y = Rel_Abundance, fill = PrimerSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "Top 10 Prey Families: Relative Abundance Comparison",
       x = "Family", y = "Relative Abundance (%)", fill = "Primer Set") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_base, "figures/03_family_abundance_comparison.pdf"), p3, width = 10, height = 6)
cat("Saved: figures/03_family_abundance_comparison.pdf\n")

# Plot: Order-only assignments
p4 <- order_only_combined %>%
  filter(Total_Reads > 0) %>%
  ggplot(aes(x = reorder(Order, -Total_Reads), y = Rel_Abundance, fill = PrimerSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "Order-Only Assignments (No Family Assignment)",
       subtitle = "Potential missed family-level assignments",
       x = "Order", y = "Relative Abundance (%)", fill = "Primer Set") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_base, "figures/03_order_only_breakdown.pdf"), p4, width = 10, height = 6)
cat("Saved: figures/03_order_only_breakdown.pdf\n")

# ==============================================================================
# Analysis 4: Poorly-Assigned ASV Deep Dive
# ==============================================================================

cat("\n=== Analysis 4: Poorly-Assigned ASV Deep Dive ===\n")

identify_poorly_assigned <- function(tax_df, table_df, seqs, primer_set, top_n = 20) {
  # Calculate total abundance per ASV
  asv_abundance <- table_df %>%
    mutate(Total_Reads = rowSums(select(., -ASV))) %>%
    select(ASV, Total_Reads)

  # Merge with taxonomy
  merged <- tax_df %>%
    inner_join(asv_abundance, by = c("Feature.ID" = "ASV"))

  # Identify poorly assigned
  poorly_assigned <- merged %>%
    filter(is.na(Family) | Confidence < 0.85) %>%
    arrange(desc(Total_Reads)) %>%
    slice_head(n = top_n) %>%
    mutate(PrimerSet = primer_set)

  # Extract sequences
  seq_subset <- seqs[names(seqs) %in% poorly_assigned$Feature.ID]

  return(list(data = poorly_assigned, seqs = seq_subset))
}

ps1_poorly <- identify_poorly_assigned(ps1_tax, ps1_table, ps1_seqs, "PS1")
ps2_poorly <- identify_poorly_assigned(ps2_tax, ps2_table, ps2_seqs, "PS2")

# Save tables
write_csv(ps1_poorly$data, file.path(output_base, "tables/04_poorly_assigned_ps1_top20.csv"))
write_csv(ps2_poorly$data, file.path(output_base, "tables/04_poorly_assigned_ps2_top20.csv"))
cat("Saved: tables/04_poorly_assigned_ps1_top20.csv\n")
cat("Saved: tables/04_poorly_assigned_ps2_top20.csv\n")

# Save sequences as FASTA
write_fasta <- function(seqs, filename) {
  fasta_lines <- character()
  for (i in seq_along(seqs)) {
    fasta_lines <- c(fasta_lines, paste0(">", names(seqs)[i]), as.character(seqs[i]))
  }
  writeLines(fasta_lines, filename)
}

write_fasta(ps1_poorly$seqs, file.path(output_base, "sequences/ps1_poorly_assigned_top20.fasta"))
write_fasta(ps2_poorly$seqs, file.path(output_base, "sequences/ps2_poorly_assigned_top20.fasta"))
cat("Saved: sequences/ps1_poorly_assigned_top20.fasta\n")
cat("Saved: sequences/ps2_poorly_assigned_top20.fasta\n")

# Summary: % of total reads that are poorly assigned
ps1_total_reads <- sum(ps1_table %>% select(-ASV))
ps2_total_reads <- sum(ps2_table %>% select(-ASV))
ps1_poorly_reads <- sum(ps1_poorly$data$Total_Reads)
ps2_poorly_reads <- sum(ps2_poorly$data$Total_Reads)

poorly_summary <- tibble(
  PrimerSet = c("PS1", "PS2"),
  Total_ASVs = c(nrow(ps1_tax), nrow(ps2_tax)),
  Poorly_Assigned_ASVs = c(
    sum(is.na(ps1_tax$Family) | ps1_tax$Confidence < 0.85),
    sum(is.na(ps2_tax$Family) | ps2_tax$Confidence < 0.85)
  ),
  Percent_ASVs_Poorly_Assigned = 100 * Poorly_Assigned_ASVs / Total_ASVs,
  Top20_Reads = c(ps1_poorly_reads, ps2_poorly_reads),
  Total_Reads = c(ps1_total_reads, ps2_total_reads),
  Percent_Reads_Top20 = 100 * Top20_Reads / Total_Reads
)
write_csv(poorly_summary, file.path(output_base, "tables/04_poorly_assigned_summary.csv"))
cat("Saved: tables/04_poorly_assigned_summary.csv\n")

# Plot
p5 <- poorly_summary %>%
  select(PrimerSet, Percent_ASVs_Poorly_Assigned, Percent_Reads_Top20) %>%
  pivot_longer(-PrimerSet, names_to = "Metric", values_to = "Percent") %>%
  mutate(Metric = recode(Metric,
                         Percent_ASVs_Poorly_Assigned = "% ASVs Poorly Assigned",
                         Percent_Reads_Top20 = "% Reads in Top20 Poorly Assigned")) %>%
  ggplot(aes(x = PrimerSet, y = Percent, fill = PrimerSet)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "Poorly-Assigned ASVs: Proportion of Dataset",
       x = "Primer Set", y = "Percent (%)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(output_base, "figures/04_poorly_assigned_abundance.pdf"), p5, width = 8, height = 4)
cat("Saved: figures/04_poorly_assigned_abundance.pdf\n")

# ==============================================================================
# Analysis 5: Paired Sample Analysis
# ==============================================================================

cat("\n=== Analysis 5: Paired Sample Analysis ===\n")

# Identify paired samples
ps1_paired_samples <- ps1_metadata %>%
  filter(both_primersets == "yes") %>%
  pull(sampleid)

ps2_paired_samples <- ps2_metadata %>%
  filter(both_primersets == "yes") %>%
  pull(sampleid) %>%
  str_remove("BF2$")

cat("PS1 paired samples:", length(ps1_paired_samples), "\n")
cat("PS2 paired samples:", length(ps2_paired_samples), "\n")

# Extract paired sample data
ps1_paired_data <- per_sample_rates %>%
  filter(PrimerSet == "PS1", Sample %in% ps1_paired_samples) %>%
  mutate(Base_Sample = Sample) %>%
  select(Base_Sample, PS1_Percent = Percent_Assigned, PS1_Total_Reads = Total_Reads)

ps2_paired_data <- per_sample_rates %>%
  filter(PrimerSet == "PS2") %>%
  mutate(Base_Sample = str_remove(Sample, "BF2$")) %>%
  filter(Base_Sample %in% ps1_paired_samples) %>%
  select(Base_Sample, PS2_Percent = Percent_Assigned, PS2_Total_Reads = Total_Reads)

paired_comparison <- ps1_paired_data %>%
  inner_join(ps2_paired_data, by = "Base_Sample") %>%
  mutate(
    Difference = PS2_Percent - PS1_Percent,
    Outlier = abs(Difference) > 20
  )

write_csv(paired_comparison, file.path(output_base, "tables/05_paired_sample_assignments.csv"))
cat("Saved: tables/05_paired_sample_assignments.csv\n")
cat("Paired samples analyzed:", nrow(paired_comparison), "\n")
cat("Outlier samples (>20% difference):", sum(paired_comparison$Outlier), "\n")

# Statistical tests
if (nrow(paired_comparison) > 2) {
  cor_test <- cor.test(paired_comparison$PS1_Percent, paired_comparison$PS2_Percent)
  wilcox_paired <- wilcox.test(paired_comparison$PS1_Percent, paired_comparison$PS2_Percent, paired = TRUE)

  stats_tests$paired_correlation <- tibble(
    Test = "Pearson Correlation: Paired Samples (PS1 vs PS2)",
    Statistic = cor_test$estimate,
    P_value = cor_test$p.value,
    Effect_Size = cor_test$estimate
  )

  stats_tests$paired_wilcox <- tibble(
    Test = "Wilcoxon Signed-Rank: Paired Samples (PS1 vs PS2)",
    Statistic = wilcox_paired$statistic,
    P_value = wilcox_paired$p.value,
    Effect_Size = NA_real_
  )

  # Plot
  p6 <- ggplot(paired_comparison, aes(x = PS1_Percent, y = PS2_Percent)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Outlier), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    labs(title = "Paired Sample Family Assignment Rates",
         subtitle = sprintf("r = %.3f, p = %.3f | n = %d samples",
                           cor_test$estimate, cor_test$p.value, nrow(paired_comparison)),
         x = "PS1: % Reads Assigned to Family",
         y = "PS2: % Reads Assigned to Family",
         color = "Outlier\n(>20% diff)") +
    theme_bw() +
    coord_fixed()

  ggsave(file.path(output_base, "figures/05_paired_sample_correlation.pdf"), p6, width = 7, height = 7)
  cat("Saved: figures/05_paired_sample_correlation.pdf\n")
} else {
  cat("Warning: Not enough paired samples for correlation analysis\n")
}

# ==============================================================================
# Analysis 6: Compositional Comparison
# ==============================================================================

cat("\n=== Analysis 6: Compositional Comparison ===\n")

# Aggregate family-level composition across all samples
calc_composition <- function(tax_df, table_df, primer_set) {
  merged <- tax_df %>%
    select(Feature.ID, Family) %>%
    inner_join(table_df, by = c("Feature.ID" = "ASV"))

  composition <- merged %>%
    mutate(Total_Reads = rowSums(select(., -Feature.ID, -Family))) %>%
    group_by(Family) %>%
    summarise(Total_Reads = sum(Total_Reads), .groups = "drop") %>%
    mutate(
      PrimerSet = primer_set,
      Rel_Abundance = 100 * Total_Reads / sum(Total_Reads)
    ) %>%
    arrange(desc(Total_Reads))

  return(composition)
}

ps1_composition <- calc_composition(ps1_tax, ps1_table, "PS1")
ps2_composition <- calc_composition(ps2_tax, ps2_table, "PS2")
composition_combined <- bind_rows(ps1_composition, ps2_composition)

write_csv(composition_combined, file.path(output_base, "tables/06_family_composition_comparison.csv"))
cat("Saved: tables/06_family_composition_comparison.csv\n")

# Chi-square test for composition differences
top_families_for_test <- composition_combined %>%
  group_by(Family) %>%
  summarise(Max_Reads = max(Total_Reads), .groups = "drop") %>%
  filter(!is.na(Family)) %>%
  slice_max(Max_Reads, n = 10) %>%
  pull(Family)

comp_test_data <- composition_combined %>%
  filter(Family %in% top_families_for_test) %>%
  select(Family, PrimerSet, Total_Reads) %>%
  pivot_wider(names_from = PrimerSet, values_from = Total_Reads, values_fill = 0)

comp_chi <- chisq.test(select(comp_test_data, -Family))
stats_tests$composition <- tibble(
  Test = "Chi-square: Family Composition (PS1 vs PS2)",
  Statistic = comp_chi$statistic,
  P_value = comp_chi$p.value,
  Effect_Size = sqrt(comp_chi$statistic / sum(select(comp_test_data, -Family)))
)

# Plot: Stacked bar composition
p7 <- composition_combined %>%
  filter(!is.na(Family)) %>%
  group_by(PrimerSet) %>%
  slice_max(Total_Reads, n = 15) %>%
  ungroup() %>%
  ggplot(aes(x = PrimerSet, y = Rel_Abundance, fill = reorder(Family, -Rel_Abundance))) +
  geom_bar(stat = "identity") +
  labs(title = "Family-Level Composition: Top 15 Families",
       x = "Primer Set", y = "Relative Abundance (%)", fill = "Family") +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(output_base, "figures/06_composition_stacked_bar.pdf"), p7, width = 10, height = 6)
cat("Saved: figures/06_composition_stacked_bar.pdf\n")

# Plot: Hydroptilidae vs Chironomidae direct comparison
key_families <- c("Hydroptilidae", "Chironomidae")
p8 <- composition_combined %>%
  filter(Family %in% key_families) %>%
  ggplot(aes(x = Family, y = Rel_Abundance, fill = PrimerSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", Rel_Abundance)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("PS1" = "#E69F00", "PS2" = "#56B4E9")) +
  labs(title = "Hydroptilidae vs Chironomidae: Direct Comparison",
       subtitle = "Key prey families showing primer set differences",
       x = "Family", y = "Relative Abundance (%)", fill = "Primer Set") +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave(file.path(output_base, "figures/06_hydroptilidae_vs_chironomidae.pdf"), p8, width = 8, height = 6)
cat("Saved: figures/06_hydroptilidae_vs_chironomidae.pdf\n")

# ==============================================================================
# Compile Statistical Tests Summary
# ==============================================================================

cat("\n=== Compiling Statistical Tests Summary ===\n")

stats_summary <- bind_rows(stats_tests)
write_csv(stats_summary, file.path(output_base, "tables/statistical_tests_summary.csv"))
cat("Saved: tables/statistical_tests_summary.csv\n")

# ==============================================================================
# Generate Markdown Report
# ==============================================================================

cat("\n=== Generating Markdown Report ===\n")

report_lines <- c(
  "# Taxonomy Assignment Rate Analysis: PS1 vs PS2",
  "",
  sprintf("**Analysis Date:** %s", Sys.Date()),
  sprintf("**Output Directory:** `%s`", output_base),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  sprintf("This analysis investigated whether differential taxonomy assignment rates between two COI primer sets (PS1: mlCOIintF+BR2, PS2: BF2+BR1) could explain observed differences in prey composition. PS1 recovered %d ASVs, while PS2 recovered %d ASVs.", nrow(ps1_tax), nrow(ps2_tax)),
  "",
  "### Key Findings",
  "",
  sprintf("- **ASV-level Family assignment:** PS1: %.1f%%, PS2: %.1f%% (p = %.3f)",
          filter(ps1_rates, Rank == "Family")$Percent_Assigned,
          filter(ps2_rates, Rank == "Family")$Percent_Assigned,
          stats_summary$P_value[stats_summary$Test == "Chi-square: ASV Family Assignment (PS1 vs PS2)"]),
  sprintf("- **Read-weighted Family assignment:** PS1: %.1f%%, PS2: %.1f%% (p = %.3f)",
          ps1_read_rates$global, ps2_read_rates$global,
          stats_summary$P_value[stats_summary$Test == "Wilcoxon: Per-Sample Family Assignment (PS1 vs PS2)"]),
  sprintf("- **Hydroptilidae abundance:** PS1: %.1f%%, PS2: %.1f%%",
          filter(ps1_composition, Family == "Hydroptilidae")$Rel_Abundance[1],
          filter(ps2_composition, Family == "Hydroptilidae")$Rel_Abundance[1]),
  sprintf("- **Chironomidae abundance:** PS1: %.1f%%, PS2: %.1f%%",
          filter(ps1_composition, Family == "Chironomidae")$Rel_Abundance[1],
          filter(ps2_composition, Family == "Chironomidae")$Rel_Abundance[1]),
  "",
  "---",
  "",
  "## 1. Global Assignment Rates",
  "",
  "### ASV-Level Assignment",
  "",
  "Assignment rates by taxonomic rank (proportion of ASVs assigned):",
  "",
  "| Rank | PS1 | PS2 |",
  "|------|-----|-----|",
  paste0("| Family | ", sprintf("%.1f%%", filter(ps1_rates, Rank == "Family")$Percent_Assigned),
         " | ", sprintf("%.1f%%", filter(ps2_rates, Rank == "Family")$Percent_Assigned), " |"),
  paste0("| Genus | ", sprintf("%.1f%%", filter(ps1_rates, Rank == "Genus")$Percent_Assigned),
         " | ", sprintf("%.1f%%", filter(ps2_rates, Rank == "Genus")$Percent_Assigned), " |"),
  paste0("| Species | ", sprintf("%.1f%%", filter(ps1_rates, Rank == "Species")$Percent_Assigned),
         " | ", sprintf("%.1f%%", filter(ps2_rates, Rank == "Species")$Percent_Assigned), " |"),
  "",
  "**See:** `figures/01_assignment_rates_by_rank.pdf`",
  "",
  "### Read-Weighted Assignment",
  "",
  sprintf("When weighted by read abundance, PS1 assigned %.1f%% of reads to Family level, while PS2 assigned %.1f%%.", ps1_read_rates$global, ps2_read_rates$global),
  "",
  "**See:** `figures/02_read_weighted_boxplot.pdf`",
  "",
  "---",
  "",
  "## 2. Taxon-Specific Analysis",
  "",
  "### Top Prey Families",
  "",
  "The 5 most abundant families by total read count:",
  "",
  "**PS1:**"
)

# Add top 5 families for PS1
top5_ps1 <- ps1_composition %>% slice_head(n = 5)
for (i in 1:nrow(top5_ps1)) {
  report_lines <- c(report_lines,
                   sprintf("%d. %s: %.1f%% (%d reads)",
                          i,
                          ifelse(is.na(top5_ps1$Family[i]), "Unassigned", top5_ps1$Family[i]),
                          top5_ps1$Rel_Abundance[i],
                          top5_ps1$Total_Reads[i]))
}

report_lines <- c(report_lines, "", "**PS2:**")

# Add top 5 families for PS2
top5_ps2 <- ps2_composition %>% slice_head(n = 5)
for (i in 1:nrow(top5_ps2)) {
  report_lines <- c(report_lines,
                   sprintf("%d. %s: %.1f%% (%d reads)",
                          i,
                          ifelse(is.na(top5_ps2$Family[i]), "Unassigned", top5_ps2$Family[i]),
                          top5_ps2$Rel_Abundance[i],
                          top5_ps2$Total_Reads[i]))
}

report_lines <- c(report_lines,
  "",
  "**See:** `figures/03_family_abundance_comparison.pdf`",
  "",
  "### Order-Only Assignments",
  "",
  "ASVs assigned to Order but not Family may represent missed taxonomic assignments:",
  "",
  sprintf("- **PS1:** %.1f%% of reads in %d order-only ASVs",
          sum(filter(order_only_combined, PrimerSet == "PS1")$Rel_Abundance),
          sum(filter(order_only_combined, PrimerSet == "PS1")$N_ASVs)),
  sprintf("- **PS2:** %.1f%% of reads in %d order-only ASVs",
          sum(filter(order_only_combined, PrimerSet == "PS2")$Rel_Abundance),
          sum(filter(order_only_combined, PrimerSet == "PS2")$N_ASVs)),
  "",
  "**See:** `figures/03_order_only_breakdown.pdf`, `tables/03_order_only_assignments.csv`",
  "",
  "---",
  "",
  "## 3. Poorly-Assigned ASVs",
  "",
  sprintf("Poorly-assigned ASVs (no Family assignment OR confidence < 0.85) represent %.1f%% of PS1 ASVs and %.1f%% of PS2 ASVs.",
          filter(poorly_summary, PrimerSet == "PS1")$Percent_ASVs_Poorly_Assigned,
          filter(poorly_summary, PrimerSet == "PS2")$Percent_ASVs_Poorly_Assigned),
  "",
  sprintf("The top 20 most abundant poorly-assigned ASVs account for %.1f%% of PS1 reads and %.1f%% of PS2 reads.",
          filter(poorly_summary, PrimerSet == "PS1")$Percent_Reads_Top20,
          filter(poorly_summary, PrimerSet == "PS2")$Percent_Reads_Top20),
  "",
  "**For manual BLAST verification:**",
  "- `sequences/ps1_poorly_assigned_top20.fasta`",
  "- `sequences/ps2_poorly_assigned_top20.fasta`",
  "",
  "**See:** `figures/04_poorly_assigned_abundance.pdf`, `tables/04_poorly_assigned_summary.csv`",
  "",
  "---",
  "",
  "## 4. Paired Sample Analysis"
)

if (nrow(paired_comparison) > 2) {
  cor_val <- stats_summary$Statistic[stats_summary$Test == "Pearson Correlation: Paired Samples (PS1 vs PS2)"]
  cor_p <- stats_summary$P_value[stats_summary$Test == "Pearson Correlation: Paired Samples (PS1 vs PS2)"]

  report_lines <- c(report_lines,
    "",
    sprintf("Analyzed %d samples sequenced with both primer sets.", nrow(paired_comparison)),
    "",
    sprintf("- **Correlation:** r = %.3f, p = %.3f", cor_val, cor_p),
    sprintf("- **Outlier samples (>20%% difference):** %d", sum(paired_comparison$Outlier)),
    "",
    "**See:** `figures/05_paired_sample_correlation.pdf`, `tables/05_paired_sample_assignments.csv`"
  )
} else {
  report_lines <- c(report_lines,
    "",
    "Insufficient paired samples for correlation analysis.",
    ""
  )
}

report_lines <- c(report_lines,
  "",
  "---",
  "",
  "## 5. Compositional Comparison",
  "",
  "### Hydroptilidae vs Chironomidae",
  ""
)

hydro_ps1 <- filter(ps1_composition, Family == "Hydroptilidae")$Rel_Abundance[1]
hydro_ps2 <- filter(ps2_composition, Family == "Hydroptilidae")$Rel_Abundance[1]
chiro_ps1 <- filter(ps1_composition, Family == "Chironomidae")$Rel_Abundance[1]
chiro_ps2 <- filter(ps2_composition, Family == "Chironomidae")$Rel_Abundance[1]

report_lines <- c(report_lines,
  "| Family | PS1 | PS2 | Fold Difference |",
  "|--------|-----|-----|-----------------|",
  sprintf("| Hydroptilidae | %.1f%% | %.1f%% | %.1fx |", hydro_ps1, hydro_ps2, hydro_ps1/hydro_ps2),
  sprintf("| Chironomidae | %.1f%% | %.1f%% | %.1fx |", chiro_ps1, chiro_ps2, chiro_ps2/chiro_ps1),
  "",
  "**See:** `figures/06_hydroptilidae_vs_chironomidae.pdf`, `figures/06_composition_stacked_bar.pdf`",
  "",
  "---",
  "",
  "## 6. Conclusions and Recommendations",
  "",
  "### Main Conclusion",
  ""
)

# Determine conclusion based on assignment rates
if (abs(ps1_read_rates$global - ps2_read_rates$global) < 5) {
  conclusion <- "Assignment rates are similar between primer sets (difference < 5%). The observed compositional differences (Hydroptilidae-dominated in PS1 vs Chironomidae-dominated in PS2) are likely REAL BIOLOGICAL DIFFERENCES in primer amplification specificity, not taxonomic assignment artifacts."
} else if (ps1_read_rates$global < ps2_read_rates$global - 10) {
  conclusion <- "PS1 has substantially lower assignment rates than PS2. The Hydroptilidae-dominated pattern in PS1 may be partially driven by TAXONOMIC ASSIGNMENT FAILURE in other groups (especially Chironomidae). Manual BLAST verification of poorly-assigned ASVs is recommended."
} else {
  conclusion <- "Assignment rates show moderate differences. The compositional differences likely reflect a COMBINATION of real biological differences in primer specificity and some taxonomic assignment variation."
}

report_lines <- c(report_lines,
  conclusion,
  "",
  "### Next Steps",
  "",
  "1. **Manual BLAST verification:** Submit poorly-assigned FASTA files to NCBI BLAST to identify potential missed assignments",
  "2. **Primer specificity analysis:** Investigate in silico primer binding sites for Hydroptilidae vs Chironomidae",
  "3. **Reference database evaluation:** Check COInr coverage for key prey families in this geographic region",
  "4. **Statistical modeling:** Use compositional data analysis (e.g., ALDEx2, ANCOM) to formally test for differential abundance",
  "",
  "---",
  "",
  "## Output Files",
  "",
  "### Tables",
  "- `01_asv_assignment_rates.csv` - ASV-level assignment by rank",
  "- `02_read_weighted_assignment_global.csv` - Global read-weighted rates",
  "- `02_read_weighted_assignment_per_sample.csv` - Per-sample rates",
  "- `03_family_specific_stats.csv` - Family-level abundances",
  "- `03_order_only_assignments.csv` - Order-only ASVs",
  "- `04_poorly_assigned_ps1_top20.csv` - Top 20 poorly-assigned PS1",
  "- `04_poorly_assigned_ps2_top20.csv` - Top 20 poorly-assigned PS2",
  "- `04_poorly_assigned_summary.csv` - Summary statistics",
  "- `05_paired_sample_assignments.csv` - Paired sample comparison",
  "- `06_family_composition_comparison.csv` - Family composition",
  "- `statistical_tests_summary.csv` - All statistical tests",
  "",
  "### Sequences",
  "- `ps1_poorly_assigned_top20.fasta` - FASTA for BLAST",
  "- `ps2_poorly_assigned_top20.fasta` - FASTA for BLAST",
  "",
  "### Figures",
  "- `01_assignment_rates_by_rank.pdf`",
  "- `02_read_weighted_boxplot.pdf`",
  "- `03_family_abundance_comparison.pdf`",
  "- `03_order_only_breakdown.pdf`",
  "- `04_poorly_assigned_abundance.pdf`",
  "- `05_paired_sample_correlation.pdf`",
  "- `06_composition_stacked_bar.pdf`",
  "- `06_hydroptilidae_vs_chironomidae.pdf`",
  "",
  "---",
  "",
  sprintf("**Report generated:** %s", Sys.time())
)

writeLines(report_lines, file.path(output_base, "assignment_diagnostics_report.md"))
cat("Saved: assignment_diagnostics_report.md\n")

cat("\n=== Analysis Complete ===\n")
cat("All outputs saved to:", output_base, "\n")
cat("Review the markdown report for detailed findings.\n")
