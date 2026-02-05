#!/usr/bin/env Rscript
# filtering_loss_analysis.R
#
# Quantify taxonomy filtering loss in gut content pipeline
# Determines whether differential filtering rates between primer sets
# could explain observed compositional differences
#
# Author: Claude Code
# Date: 2026-02-04

# ============================================================================
# SETUP AND CONFIGURATION
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(qiime2R)
  library(patchwork)
})

# Get root directory from environment
root <- Sys.getenv("root")
if (root == "") {
  stop("ERROR: 'root' environment variable not set. Source gut_contents.env first.")
}

# Define output directory
output_base <- file.path(root, "output_data", "07_Analysis", "filtering_loss_analysis")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "report"), recursive = TRUE, showWarnings = FALSE)

# Define input file paths
input_files <- list(
  ps1_before = file.path(root, "output_data", "03_Clustered_Data", "primerset1_all_p985_table.qza"),
  ps1_after = file.path(root, "output_data", "06_Generate_Output", "primerset1_all_p985_table_filtd.qza"),
  ps1_taxonomy = file.path(root, "output_data", "06_Generate_Output", "primerset1_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv"),
  ps2_before = file.path(root, "output_data", "03_Clustered_Data", "primerset2_all_p985_table.qza"),
  ps2_after = file.path(root, "output_data", "06_Generate_Output", "primerset2_all_p985_table_filtd.qza"),
  ps2_taxonomy = file.path(root, "output_data", "06_Generate_Output", "primerset2_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv")
)

# Verify all input files exist
for (file_name in names(input_files)) {
  if (!file.exists(input_files[[file_name]])) {
    stop(paste0("ERROR: Input file not found: ", input_files[[file_name]]))
  }
}

cat("✓ Output directory created:", output_base, "\n")
cat("✓ All input files verified\n\n")

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

#' Load QIIME2 feature table from .qza file
#'
#' @param qza_file Path to .qza file
#' @return Tibble with ASV column and sample columns
load_feature_table_qza <- function(qza_file) {
  cat("  Loading:", basename(qza_file), "... ")

  # Read QZA file
  qza_data <- qiime2R::read_qza(qza_file)

  # Extract feature table (stored in $data)
  feature_table <- qza_data$data

  # Convert to tibble with ASV as first column
  feature_tibble <- feature_table %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    as_tibble()

  # Calculate total reads
  total_reads <- sum(feature_table)
  n_asvs <- nrow(feature_tibble)
  n_samples <- ncol(feature_tibble) - 1

  cat(sprintf("✓ (%d ASVs, %d samples, %s reads)\n",
              n_asvs, n_samples, format(total_reads, big.mark = ",")))

  return(feature_tibble)
}

#' Parse VSEARCH taxonomy string
#'
#' @param taxon_string Taxonomy string from VSEARCH
#' @return Named list with Kingdom, Phylum, Class, Order, Family, Genus, Species
parse_vsearch_taxonomy <- function(taxon_string) {
  # Handle "Unassigned" case
  if (taxon_string == "Unassigned" || is.na(taxon_string)) {
    return(list(Kingdom = NA_character_, Phylum = NA_character_,
                Class = NA_character_, Order = NA_character_,
                Family = NA_character_, Genus = NA_character_,
                Species = NA_character_))
  }

  # Split on "; "
  ranks <- str_split(taxon_string, "; ")[[1]]

  # Initialize result
  result <- list(Kingdom = NA_character_, Phylum = NA_character_,
                 Class = NA_character_, Order = NA_character_,
                 Family = NA_character_, Genus = NA_character_,
                 Species = NA_character_)

  # Rank prefix mapping
  rank_map <- c("k__" = "Kingdom", "p__" = "Phylum", "c__" = "Class",
                "o__" = "Order", "f__" = "Family", "g__" = "Genus", "s__" = "Species")

  for (rank_str in ranks) {
    # Extract prefix
    prefix <- str_extract(rank_str, "^[kpcofgs]__")

    if (!is.na(prefix) && prefix %in% names(rank_map)) {
      rank_name <- rank_map[prefix]

      # Remove prefix and taxid suffix
      value <- str_remove(rank_str, "^[kpcofgs]__")
      value <- str_remove(value, "_[0-9]+$")
      value <- str_trim(value)

      # Convert empty to NA
      if (value != "" && !str_detect(value, "^[kpcofgs]__$")) {
        result[[rank_name]] <- value
      }
    }
  }

  return(result)
}

#' Load and parse raw taxonomy TSV file
#'
#' @param tsv_file Path to taxonomy TSV file
#' @return Tibble with ASV and taxonomic ranks
load_raw_taxonomy_tsv <- function(tsv_file) {
  cat("  Loading:", basename(tsv_file), "... ")

  # Read TSV file
  taxonomy_raw <- read_tsv(tsv_file, col_types = cols(.default = "c"))

  # Check required columns
  if (!all(c("Feature ID", "Taxon") %in% colnames(taxonomy_raw))) {
    stop("ERROR: Taxonomy file missing required columns")
  }

  # Parse taxonomy strings
  taxonomy_parsed <- taxonomy_raw %>%
    rename(ASV = `Feature ID`) %>%
    rowwise() %>%
    mutate(
      parsed = list(parse_vsearch_taxonomy(Taxon))
    ) %>%
    ungroup() %>%
    mutate(
      Kingdom = map_chr(parsed, ~ .x$Kingdom %||% NA_character_),
      Phylum = map_chr(parsed, ~ .x$Phylum %||% NA_character_),
      Class = map_chr(parsed, ~ .x$Class %||% NA_character_),
      Order = map_chr(parsed, ~ .x$Order %||% NA_character_),
      Family = map_chr(parsed, ~ .x$Family %||% NA_character_),
      Genus = map_chr(parsed, ~ .x$Genus %||% NA_character_),
      Species = map_chr(parsed, ~ .x$Species %||% NA_character_)
    ) %>%
    select(-parsed)

  n_asvs <- nrow(taxonomy_parsed)
  cat(sprintf("✓ (%d ASVs)\n", n_asvs))

  return(taxonomy_parsed)
}

# ============================================================================
# LOAD ALL DATA
# ============================================================================

cat("Loading data files...\n")

# Load PS1 data
cat("\nPrimer Set 1:\n")
ps1_before <- load_feature_table_qza(input_files$ps1_before)
ps1_after <- load_feature_table_qza(input_files$ps1_after)
ps1_taxonomy <- load_raw_taxonomy_tsv(input_files$ps1_taxonomy)

# Load PS2 data
cat("\nPrimer Set 2:\n")
ps2_before <- load_feature_table_qza(input_files$ps2_before)
ps2_after <- load_feature_table_qza(input_files$ps2_after)
ps2_taxonomy <- load_raw_taxonomy_tsv(input_files$ps2_taxonomy)

# Validate that after is subset of before
cat("\nValidating data integrity...\n")
ps1_after_in_before <- all(ps1_after$ASV %in% ps1_before$ASV)
ps2_after_in_before <- all(ps2_after$ASV %in% ps2_before$ASV)

if (!ps1_after_in_before) {
  stop("ERROR: PS1 after table contains ASVs not in before table")
}
if (!ps2_after_in_before) {
  stop("ERROR: PS2 after table contains ASVs not in before table")
}

cat("✓ After tables are proper subsets of before tables\n")

# ============================================================================
# ANALYSIS 1: GLOBAL FILTERING LOSS
# ============================================================================

#' Calculate global filtering loss statistics
#'
#' @param before_table Feature table before filtering
#' @param after_table Feature table after filtering
#' @param primerset Primer set name
#' @return Tibble with summary statistics
calc_global_filtering_loss <- function(before_table, after_table, primerset) {
  # Count ASVs
  n_asvs_before <- nrow(before_table)
  n_asvs_after <- nrow(after_table)
  n_asvs_lost <- n_asvs_before - n_asvs_after
  pct_asvs_lost <- 100 * n_asvs_lost / n_asvs_before

  # Sum reads
  total_reads_before <- sum(select(before_table, -ASV))
  total_reads_after <- sum(select(after_table, -ASV))
  total_reads_lost <- total_reads_before - total_reads_after
  pct_reads_lost <- 100 * total_reads_lost / total_reads_before

  tibble(
    PrimerSet = primerset,
    ASVs_Before = n_asvs_before,
    ASVs_After = n_asvs_after,
    ASVs_Lost = n_asvs_lost,
    Percent_ASVs_Lost = pct_asvs_lost,
    Reads_Before = total_reads_before,
    Reads_After = total_reads_after,
    Reads_Lost = total_reads_lost,
    Percent_Reads_Lost = pct_reads_lost
  )
}

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS 1: Global Filtering Loss\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

global_summary <- bind_rows(
  calc_global_filtering_loss(ps1_before, ps1_after, "PS1"),
  calc_global_filtering_loss(ps2_before, ps2_after, "PS2")
)

print(global_summary)

# Save table
write_csv(global_summary, file.path(output_base, "tables", "global_filtering_summary.csv"))
cat("\n✓ Saved: tables/global_filtering_summary.csv\n")

# ============================================================================
# ANALYSIS 2: PER-SAMPLE FILTERING LOSS
# ============================================================================

#' Calculate per-sample filtering loss
#'
#' @param before_table Feature table before filtering
#' @param after_table Feature table after filtering
#' @param primerset Primer set name
#' @return Long-format tibble with per-sample statistics
calc_per_sample_loss <- function(before_table, after_table, primerset) {
  # Get sample columns (all except ASV)
  sample_cols <- setdiff(colnames(before_table), "ASV")

  # Calculate reads per sample before and after
  per_sample_stats <- tibble(
    Sample = sample_cols,
    PrimerSet = primerset,
    Reads_Before = map_dbl(sample_cols, ~ sum(before_table[[.x]])),
    Reads_After = map_dbl(sample_cols, ~ sum(after_table[[.x]]))
  ) %>%
    mutate(
      Reads_Lost = Reads_Before - Reads_After,
      Percent_Lost = 100 * Reads_Lost / Reads_Before
    )

  return(per_sample_stats)
}

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS 2: Per-Sample Filtering Loss\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

per_sample_summary <- bind_rows(
  calc_per_sample_loss(ps1_before, ps1_after, "PS1"),
  calc_per_sample_loss(ps2_before, ps2_after, "PS2")
)

# Display summary statistics
per_sample_stats <- per_sample_summary %>%
  group_by(PrimerSet) %>%
  summarise(
    N_Samples = n(),
    Mean_Percent_Lost = mean(Percent_Lost),
    Median_Percent_Lost = median(Percent_Lost),
    SD_Percent_Lost = sd(Percent_Lost),
    Min_Percent_Lost = min(Percent_Lost),
    Max_Percent_Lost = max(Percent_Lost),
    .groups = "drop"
  )

print(per_sample_stats)

# Save table
write_csv(per_sample_summary, file.path(output_base, "tables", "per_sample_filtering_loss.csv"))
cat("\n✓ Saved: tables/per_sample_filtering_loss.csv\n")

# Warning for high-loss samples
high_loss_samples <- per_sample_summary %>%
  filter(Percent_Lost > 50)

if (nrow(high_loss_samples) > 0) {
  cat("\n⚠ WARNING:", nrow(high_loss_samples), "samples lost >50% of reads during filtering\n")
}

# ============================================================================
# ANALYSIS 3: CATEGORIZE REMOVED ASVs
# ============================================================================

#' Categorize removed ASV by taxonomic assignment level
#'
#' @param Kingdom Kingdom assignment
#' @param Phylum Phylum assignment
#' @param Class Class assignment
#' @param Order Order assignment
#' @param Family Family assignment
#' @return Category string
categorize_removed_asv <- function(Kingdom, Phylum, Class, Order, Family) {
  # Decision tree for categorization
  if (is.na(Kingdom) && is.na(Phylum) && is.na(Class) && is.na(Order) && is.na(Family)) {
    return("Unassigned")
  } else if (!is.na(Order) && is.na(Family)) {
    return("Order-only (no Family)")
  } else if (!is.na(Class) && is.na(Order)) {
    return("Class-only (no Order)")
  } else if (!is.na(Phylum) && is.na(Class)) {
    return("Phylum-only")
  } else if (!is.na(Kingdom) && is.na(Phylum)) {
    return("Kingdom-only")
  } else {
    return("Other")
  }
}

#' Characterize removed ASVs
#'
#' @param before_table Feature table before filtering
#' @param after_table Feature table after filtering
#' @param taxonomy Taxonomy table
#' @param primerset Primer set name
#' @return Tibble with removed ASV details
characterize_removed_asvs <- function(before_table, after_table, taxonomy, primerset) {
  # Identify removed ASVs
  removed_asv_ids <- setdiff(before_table$ASV, after_table$ASV)

  # Get total reads for each removed ASV
  removed_asv_reads <- before_table %>%
    filter(ASV %in% removed_asv_ids) %>%
    rowwise() %>%
    mutate(Total_Reads = sum(c_across(-ASV))) %>%
    ungroup() %>%
    select(ASV, Total_Reads)

  # Join with taxonomy and categorize
  removed_asv_details <- removed_asv_reads %>%
    left_join(taxonomy, by = "ASV") %>%
    rowwise() %>%
    mutate(
      Category = categorize_removed_asv(Kingdom, Phylum, Class, Order, Family)
    ) %>%
    ungroup() %>%
    mutate(PrimerSet = primerset)

  return(removed_asv_details)
}

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS 3: Taxonomy Breakdown of Removed ASVs\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

# Characterize removed ASVs for both primer sets
ps1_removed <- characterize_removed_asvs(ps1_before, ps1_after, ps1_taxonomy, "PS1")
ps2_removed <- characterize_removed_asvs(ps2_before, ps2_after, ps2_taxonomy, "PS2")

all_removed <- bind_rows(ps1_removed, ps2_removed)

# Summarize by category
category_summary <- all_removed %>%
  group_by(PrimerSet, Category) %>%
  summarise(
    N_ASVs = n(),
    Total_Reads = sum(Total_Reads),
    .groups = "drop"
  ) %>%
  group_by(PrimerSet) %>%
  mutate(
    Percent_of_Lost_ASVs = 100 * N_ASVs / sum(N_ASVs),
    Percent_of_Lost_Reads = 100 * Total_Reads / sum(Total_Reads)
  ) %>%
  ungroup() %>%
  arrange(PrimerSet, desc(Total_Reads))

print(category_summary)

# Save tables
write_csv(category_summary, file.path(output_base, "tables", "filtered_asv_taxonomy_breakdown.csv"))
write_csv(all_removed, file.path(output_base, "tables", "removed_asvs_full_list.csv"))

cat("\n✓ Saved: tables/filtered_asv_taxonomy_breakdown.csv\n")
cat("✓ Saved: tables/removed_asvs_full_list.csv\n")

# ============================================================================
# ANALYSIS 4: ORDER-LEVEL LOSS
# ============================================================================

#' Analyze Order-level loss
#'
#' @param removed_asvs_df Tibble with removed ASV details
#' @return Tibble with Order-level summary
analyze_order_loss <- function(removed_asvs_df) {
  # Filter for Order-only category
  order_only <- removed_asvs_df %>%
    filter(Category == "Order-only (no Family)") %>%
    filter(!is.na(Order))

  # Summarize by Order
  order_summary <- order_only %>%
    group_by(PrimerSet, Order) %>%
    summarise(
      N_ASVs = n(),
      Total_Reads = sum(Total_Reads),
      .groups = "drop"
    ) %>%
    arrange(PrimerSet, desc(Total_Reads))

  return(order_summary)
}

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS 4: Order-Level Loss Analysis\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

order_loss_summary <- analyze_order_loss(all_removed)

# Show top 10 orders by total reads for each primer set
cat("Top 10 Orders by Lost Reads:\n\n")
cat("Primer Set 1:\n")
print(head(filter(order_loss_summary, PrimerSet == "PS1"), 10))
cat("\nPrimer Set 2:\n")
print(head(filter(order_loss_summary, PrimerSet == "PS2"), 10))

# Save table
write_csv(order_loss_summary, file.path(output_base, "tables", "lost_reads_by_order.csv"))
cat("\n✓ Saved: tables/lost_reads_by_order.csv\n")

# ============================================================================
# ANALYSIS 5: STATISTICAL TESTS
# ============================================================================

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS 5: Statistical Tests\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

# Test 1: Two-proportion z-test for global read loss
ps1_global <- filter(global_summary, PrimerSet == "PS1")
ps2_global <- filter(global_summary, PrimerSet == "PS2")

prop_test_result <- prop.test(
  x = c(ps1_global$Reads_Lost, ps2_global$Reads_Lost),
  n = c(ps1_global$Reads_Before, ps2_global$Reads_Before)
)

# Test 2: Wilcoxon rank-sum test for per-sample loss distribution
ps1_per_sample <- filter(per_sample_summary, PrimerSet == "PS1")
ps2_per_sample <- filter(per_sample_summary, PrimerSet == "PS2")

wilcox_test_result <- wilcox.test(
  ps1_per_sample$Percent_Lost,
  ps2_per_sample$Percent_Lost
)

# Test 3: Chi-square test for category distribution
category_contingency <- category_summary %>%
  select(PrimerSet, Category, N_ASVs) %>%
  pivot_wider(names_from = PrimerSet, values_from = N_ASVs, values_fill = 0) %>%
  column_to_rownames("Category") %>%
  as.matrix()

chisq_test_result <- chisq.test(category_contingency)

# Compile test results
statistical_tests <- tibble(
  Test = c(
    "Two-proportion z-test (Global % Reads Lost)",
    "Wilcoxon rank-sum test (Per-sample % Reads Lost)",
    "Chi-square test (Category Distribution)"
  ),
  Statistic = c(
    prop_test_result$statistic,
    wilcox_test_result$statistic,
    chisq_test_result$statistic
  ),
  P_value = c(
    prop_test_result$p.value,
    wilcox_test_result$p.value,
    chisq_test_result$p.value
  ),
  Interpretation = c(
    ifelse(prop_test_result$p.value < 0.05, "Significant difference", "No significant difference"),
    ifelse(wilcox_test_result$p.value < 0.05, "Significant difference", "No significant difference"),
    ifelse(chisq_test_result$p.value < 0.05, "Category distribution differs", "Category distribution similar")
  )
)

print(statistical_tests)

# Save table
write_csv(statistical_tests, file.path(output_base, "tables", "statistical_tests_summary.csv"))
cat("\n✓ Saved: tables/statistical_tests_summary.csv\n")

# Warning for large differences
ps1_loss_pct <- ps1_global$Percent_Reads_Lost
ps2_loss_pct <- ps2_global$Percent_Reads_Lost
loss_diff <- abs(ps1_loss_pct - ps2_loss_pct)

if (loss_diff > 10) {
  cat("\n⚠ WARNING: PS1 and PS2 show >10% difference in filtering loss rate\n")
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

# Define color palette
primer_colors <- c("PS1" = "#E69F00", "PS2" = "#56B4E9")

#' Plot global filtering loss comparison
plot_filtering_loss_comparison <- function(global_summary) {
  # Prepare data for plotting
  plot_data <- global_summary %>%
    select(PrimerSet, Percent_ASVs_Lost, Percent_Reads_Lost) %>%
    pivot_longer(
      cols = c(Percent_ASVs_Lost, Percent_Reads_Lost),
      names_to = "Metric",
      values_to = "Percent_Lost"
    ) %>%
    mutate(
      Metric = recode(Metric,
                     "Percent_ASVs_Lost" = "% ASVs Lost",
                     "Percent_Reads_Lost" = "% Reads Lost")
    )

  # Create plot
  p <- ggplot(plot_data, aes(x = PrimerSet, y = Percent_Lost, fill = PrimerSet)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Percent_Lost)),
              vjust = -0.5, size = 4) +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_fill_manual(values = primer_colors) +
    labs(
      title = "Global Filtering Loss Comparison",
      subtitle = "ASVs and reads lost during Family-level filtering",
      x = "Primer Set",
      y = "Percent Lost (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90")
    )

  return(p)
}

#' Plot per-sample loss distributions
plot_per_sample_distributions <- function(per_sample_data, wilcox_pvalue) {
  # Boxplot
  p1 <- ggplot(per_sample_data, aes(x = PrimerSet, y = Percent_Lost, fill = PrimerSet)) +
    geom_boxplot(width = 0.5, outlier.shape = 21) +
    scale_fill_manual(values = primer_colors) +
    labs(
      title = "Per-Sample Filtering Loss Distribution",
      subtitle = sprintf("Wilcoxon test p = %.4f", wilcox_pvalue),
      x = "Primer Set",
      y = "% Reads Lost per Sample"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )

  # Histogram with density overlay
  p2 <- ggplot(per_sample_data, aes(x = Percent_Lost, fill = PrimerSet)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values = primer_colors) +
    labs(
      x = "% Reads Lost per Sample",
      y = "Density",
      fill = "Primer Set"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top"
    )

  # Combine plots
  combined <- p1 / p2

  return(combined)
}

#' Plot taxonomy breakdown of lost ASVs
plot_taxonomy_breakdown <- function(category_summary) {
  # ASV count breakdown
  p1 <- ggplot(category_summary, aes(x = PrimerSet, y = Percent_of_Lost_ASVs, fill = Category)) +
    geom_col() +
    geom_text(
      data = category_summary %>%
        group_by(PrimerSet) %>%
        mutate(label_y = cumsum(Percent_of_Lost_ASVs) - Percent_of_Lost_ASVs / 2),
      aes(y = label_y, label = ifelse(Percent_of_Lost_ASVs > 5, sprintf("%.1f%%", Percent_of_Lost_ASVs), "")),
      size = 3
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Lost ASVs by Taxonomic Assignment Level",
      subtitle = "By ASV count",
      x = "Primer Set",
      y = "% of Lost ASVs",
      fill = "Category"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right"
    )

  # Read count breakdown
  p2 <- ggplot(category_summary, aes(x = PrimerSet, y = Percent_of_Lost_Reads, fill = Category)) +
    geom_col() +
    geom_text(
      data = category_summary %>%
        group_by(PrimerSet) %>%
        mutate(label_y = cumsum(Percent_of_Lost_Reads) - Percent_of_Lost_Reads / 2),
      aes(y = label_y, label = ifelse(Percent_of_Lost_Reads > 5, sprintf("%.1f%%", Percent_of_Lost_Reads), "")),
      size = 3
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      subtitle = "By read abundance",
      x = "Primer Set",
      y = "% of Lost Reads",
      fill = "Category"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right"
    )

  # Combine plots
  combined <- p1 + p2 + plot_layout(guides = "collect")

  return(combined)
}

#' Plot Order-level loss
plot_order_loss <- function(order_summary) {
  # Get top 10 orders by total reads across both primer sets
  top_orders <- order_summary %>%
    group_by(Order) %>%
    summarise(Total_Reads_All = sum(Total_Reads), .groups = "drop") %>%
    arrange(desc(Total_Reads_All)) %>%
    head(10) %>%
    pull(Order)

  # Filter data
  plot_data <- order_summary %>%
    filter(Order %in% top_orders) %>%
    mutate(Order = fct_reorder(Order, Total_Reads, .fun = sum))

  # Create plot
  p <- ggplot(plot_data, aes(x = Total_Reads, y = Order, fill = PrimerSet)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = primer_colors) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = "Top 10 Orders Lost to Family-Level Filtering",
      subtitle = "Orders with assignments but no Family assignment",
      x = "Total Reads Lost",
      y = "Order",
      fill = "Primer Set"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top"
    )

  return(p)
}

#' Plot abundance distribution of removed ASVs
plot_abundance_distribution <- function(removed_asvs_df) {
  # Calculate summary statistics
  summary_stats <- removed_asvs_df %>%
    group_by(PrimerSet) %>%
    summarise(
      Median = median(Total_Reads),
      Mean = mean(Total_Reads),
      .groups = "drop"
    )

  # Create plot
  p <- ggplot(removed_asvs_df, aes(x = Total_Reads, fill = PrimerSet)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
    geom_vline(
      data = summary_stats,
      aes(xintercept = Median, color = PrimerSet),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = primer_colors) +
    scale_color_manual(values = primer_colors) +
    facet_wrap(~ PrimerSet, ncol = 1) +
    labs(
      title = "Abundance Distribution of Removed ASVs",
      subtitle = "Dashed line = median abundance",
      x = "Total Reads (log scale)",
      y = "Number of ASVs",
      fill = "Primer Set"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none"
    )

  return(p)
}

# ============================================================================
# GENERATE VISUALIZATIONS
# ============================================================================

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("Generating Visualizations\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

# Figure 1: Global filtering loss
cat("  Creating: filtering_loss_comparison.pdf ... ")
fig1 <- plot_filtering_loss_comparison(global_summary)
ggsave(
  file.path(output_base, "figures", "filtering_loss_comparison.pdf"),
  fig1,
  width = 10,
  height = 5
)
cat("✓\n")

# Figure 2: Per-sample distributions
cat("  Creating: per_sample_loss_distribution.pdf ... ")
fig2 <- plot_per_sample_distributions(per_sample_summary, wilcox_test_result$p.value)
ggsave(
  file.path(output_base, "figures", "per_sample_loss_distribution.pdf"),
  fig2,
  width = 10,
  height = 10
)
cat("✓\n")

# Figure 3: Taxonomy breakdown
cat("  Creating: lost_taxonomy_breakdown.pdf ... ")
fig3 <- plot_taxonomy_breakdown(category_summary)
ggsave(
  file.path(output_base, "figures", "lost_taxonomy_breakdown.pdf"),
  fig3,
  width = 12,
  height = 6
)
cat("✓\n")

# Figure 4: Order-level loss
cat("  Creating: lost_reads_by_order.pdf ... ")
fig4 <- plot_order_loss(order_loss_summary)
ggsave(
  file.path(output_base, "figures", "lost_reads_by_order.pdf"),
  fig4,
  width = 10,
  height = 6
)
cat("✓\n")

# Figure 5: Abundance distribution
cat("  Creating: removed_asv_abundance_distribution.pdf ... ")
fig5 <- plot_abundance_distribution(all_removed)
ggsave(
  file.path(output_base, "figures", "removed_asv_abundance_distribution.pdf"),
  fig5,
  width = 10,
  height = 8
)
cat("✓\n")

# ============================================================================
# GENERATE MARKDOWN REPORT
# ============================================================================

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("Generating Markdown Report\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

# Construct report
report_lines <- c(
  "# Filtering Loss Analysis Report",
  "",
  sprintf("**Date:** %s", Sys.Date()),
  sprintf("**Pipeline Step:** Family-level taxonomic filtering (filterVsearch_nbClassifier.R)"),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  sprintf("This analysis quantifies the loss of ASVs and reads during Family-level taxonomic filtering in the gut content metabarcoding pipeline. The pipeline retains only ASVs assigned to at least **Family level** (requiring both Order AND Family to be non-NA)."),
  "",
  "### Key Findings:",
  "",
  sprintf("- **PS1** lost **%.1f%%** of ASVs (%.1f%% of reads)",
          ps1_global$Percent_ASVs_Lost, ps1_global$Percent_Reads_Lost),
  sprintf("- **PS2** lost **%.1f%%** of ASVs (%.1f%% of reads)",
          ps2_global$Percent_ASVs_Lost, ps2_global$Percent_Reads_Lost),
  sprintf("- Absolute difference: **%.1f%%** (ASVs), **%.1f%%** (reads)",
          abs(ps1_global$Percent_ASVs_Lost - ps2_global$Percent_ASVs_Lost),
          abs(ps1_global$Percent_Reads_Lost - ps2_global$Percent_Reads_Lost)),
  sprintf("- Statistical significance (global read loss): **%s** (p = %.4f)",
          ifelse(prop_test_result$p.value < 0.05, "YES", "NO"),
          prop_test_result$p.value),
  "",
  "---",
  "",
  "## 1. Global Filtering Loss",
  "",
  "### Summary Table",
  ""
)

# Add global summary table
global_table <- knitr::kable(
  global_summary %>%
    mutate(
      ASVs_Lost = sprintf("%d (%.1f%%)", ASVs_Lost, Percent_ASVs_Lost),
      Reads_Lost = sprintf("%s (%.1f%%)", format(Reads_Lost, big.mark = ","), Percent_Reads_Lost)
    ) %>%
    select(PrimerSet, ASVs_Before, ASVs_After, ASVs_Lost, Reads_Before, Reads_After, Reads_Lost),
  format = "markdown"
)

report_lines <- c(report_lines, global_table, "", "")

# Per-sample loss section
report_lines <- c(
  report_lines,
  "## 2. Per-Sample Filtering Loss",
  "",
  "### Summary Statistics",
  ""
)

per_sample_table <- knitr::kable(per_sample_stats, format = "markdown", digits = 2)
report_lines <- c(report_lines, per_sample_table, "", "")

# Taxonomy breakdown section
report_lines <- c(
  report_lines,
  "## 3. Taxonomy Breakdown of Lost ASVs",
  "",
  "### Category Summary",
  ""
)

category_table <- knitr::kable(
  category_summary %>%
    select(PrimerSet, Category, N_ASVs, Percent_of_Lost_ASVs, Total_Reads, Percent_of_Lost_Reads),
  format = "markdown",
  digits = 1
)
report_lines <- c(report_lines, category_table, "", "")

# Order-level loss section
report_lines <- c(
  report_lines,
  "## 4. Order-Level Loss Analysis",
  "",
  "### Top 10 Orders (PS1)",
  ""
)

ps1_top_orders <- order_loss_summary %>%
  filter(PrimerSet == "PS1") %>%
  head(10)

ps1_order_table <- knitr::kable(ps1_top_orders, format = "markdown", digits = 0)
report_lines <- c(report_lines, ps1_order_table, "", "")

report_lines <- c(
  report_lines,
  "### Top 10 Orders (PS2)",
  ""
)

ps2_top_orders <- order_loss_summary %>%
  filter(PrimerSet == "PS2") %>%
  head(10)

ps2_order_table <- knitr::kable(ps2_top_orders, format = "markdown", digits = 0)
report_lines <- c(report_lines, ps2_order_table, "", "")

# Statistical tests section
report_lines <- c(
  report_lines,
  "## 5. Statistical Tests",
  ""
)

stats_table <- knitr::kable(statistical_tests, format = "markdown", digits = 4)
report_lines <- c(report_lines, stats_table, "", "")

# Conclusions section
interpretation <- if (loss_diff < 5) {
  "**Interpretation:** Filtering loss rates are **similar** between primer sets (difference < 5%). Observed compositional differences are likely **real biological differences** in primer amplification specificity, not taxonomic assignment artifacts."
} else if (loss_diff < 10) {
  "**Interpretation:** Filtering loss rates show **moderate differences** between primer sets (5-10%). Some compositional differences may be partially influenced by differential filtering, but biological differences in primer specificity are likely still the primary driver."
} else {
  "**Interpretation:** Filtering loss rates show **substantial differences** between primer sets (> 10%). Compositional differences may be **partially driven by differential filtering**—some taxa may have been disproportionately removed due to poorer Family-level assignments in one primer set."
}

report_lines <- c(
  report_lines,
  "## 6. Conclusions",
  "",
  interpretation,
  "",
  "---",
  "",
  "## Output Files",
  "",
  "### Tables (CSV)",
  "- `global_filtering_summary.csv`",
  "- `per_sample_filtering_loss.csv`",
  "- `filtered_asv_taxonomy_breakdown.csv`",
  "- `lost_reads_by_order.csv`",
  "- `statistical_tests_summary.csv`",
  "- `removed_asvs_full_list.csv`",
  "",
  "### Figures (PDF)",
  "- `filtering_loss_comparison.pdf`",
  "- `per_sample_loss_distribution.pdf`",
  "- `lost_taxonomy_breakdown.pdf`",
  "- `lost_reads_by_order.pdf`",
  "- `removed_asv_abundance_distribution.pdf`",
  ""
)

# Write report
report_path <- file.path(output_base, "report", "filtering_loss_report.md")
writeLines(report_lines, report_path)

cat("✓ Saved: report/filtering_loss_report.md\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

cat("Output directory:", output_base, "\n")
cat("  - Tables: 6 CSV files\n")
cat("  - Figures: 5 PDF files\n")
cat("  - Report: filtering_loss_report.md\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Review filtering_loss_report.md for key findings\n")
cat("  2. Examine figures to visualize filtering patterns\n")
cat("  3. Check statistical_tests_summary.csv for significance tests\n")
cat("\n")
