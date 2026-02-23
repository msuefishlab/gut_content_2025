# ===== DATA LOADING =====

#' Load primer set data into phyloseq object
#'
#' @param primerset_num Integer (1 or 2) indicating which primer set
#' @param use_no_host Logical, use NO_HOST filtered data (default TRUE)
#' @return phyloseq object
load_primerset_data <- function(primerset_num, use_no_host = TRUE) {
  ps_name <- paste0("primerset", primerset_num)

  # Construct file paths
  metadata_file <- file.path(
    root, "input_data", "metadata",
    "fish_id_metadata.tsv"
  )

  if (use_no_host) {
    table_file <- file.path(
      root, "output_data", "06_Generate_Output",
      paste0(ps_name, "_all_p985_table_filtd_NO_HOST.rarefied.qza")
    )
  } else {
    table_file <- file.path(
      root, "output_data", "06_Generate_Output",
      paste0(ps_name, "_all_p985_table_filtd.rarefied.qza")
    )
  }

  taxa_file <- file.path(
    root, "output_data", "06_Generate_Output",
    paste0(ps_name, "_all_p985_taxa_filtd_ALL.qza")
  )

  # Load data
  metadata <- read_q2metadata(metadata_file)
  feature_table <- read_qza(table_file)$data
  taxonomy <- read_qza(taxa_file)$data %>% parse_taxonomy()

  # Create phyloseq object
  OTU <- otu_table(as.matrix(feature_table), taxa_are_rows = TRUE)
  TAX <- tax_table(as.matrix(taxonomy))
  META <- sample_data(metadata)
  sample_names(META) <- META$SampleID

  ps <- phyloseq(OTU, TAX, META)

  cat("Loaded", ps_name, ":", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")

  return(ps)
}


# ===== ALPHA DIVERSITY ANALYSIS =====

#' Analyze alpha diversity with appropriate statistical test
#'
#' @param ps phyloseq object
#' @param group_var Character, column name for grouping variable
#' @return List with diversity dataframe, test result, and test type
analyze_alpha_diversity <- function(ps, group_var) {
  # Calculate diversity
  div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1"))

  # Add grouping variable
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  div[[group_var]] <- meta[[group_var]]
  div[["Sample_Name"]] <- meta$Sample_Name
  # div$SampleID <- rownames(div)

  # Test normality
  shapiro_p <- shapiro.test(div$Shannon)$p.value

  # Choose appropriate test
  if (shapiro_p < 0.05) {
    # Non-parametric
    test_result <- kruskal.test(as.formula(paste("Shannon ~", group_var)),
      data = div
    )
    test_type <- "Kruskal-Wallis"
  } else {
    # Parametric ANOVA
    test_result <- aov(as.formula(paste("Shannon ~", group_var)), data = div)
    test_type <- "ANOVA"
  }

  return(list(
    diversity = div,
    test = test_result,
    test_type = test_type,
    shapiro_p = shapiro_p
  ))
}


# ===== BETA DIVERSITY ANALYSIS =====

#' Run complete beta diversity analysis (distance + PCoA + PERMANOVA)
#'
#' @param ps phyloseq object
#' @param formula_str Character, right-hand side of formula (e.g., "~ species")
#' @return List with distance matrix, PCoA results, PERMANOVA result
run_beta_diversity_analysis <- function(ps, formula_str) {
  # Extract OTU matrix (samples as rows)
  otu_mat <- t(as.matrix(otu_table(ps)))

  # Get matching metadata
  # In phyloseq, sample IDs are stored as rownames, not in a column
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

  # Ensure metadata rownames match OTU matrix rownames
  if (!identical(rownames(meta), rownames(otu_mat))) {
    meta <- meta[match(rownames(otu_mat), rownames(meta)), , drop = FALSE]
    rownames(meta) <- rownames(otu_mat)
  }

  # Calculate Bray-Curtis distance
  dist_mat <- vegdist(otu_mat, method = "bray")

  # PCoA
  pcoa_result <- cmdscale(dist_mat, eig = TRUE, k = 2)
  var_exp <- pcoa_result$eig / sum(pcoa_result$eig) * 100

  # PERMANOVA - use dist_mat directly, not in formula
  # Formula should just be the right-hand side (e.g., "~ species")
  perm_result <- adonis2(as.formula(paste("dist_mat", formula_str)),
    data = meta,
    permutations = 999
  )

  return(list(
    distance = dist_mat,
    pcoa = pcoa_result,
    variance_explained = var_exp,
    permanova = perm_result,
    metadata = meta,
    otu_matrix = otu_mat
  ))
}


# ===== COMPARATIVE PLOTTING =====

#' Plot Shannon diversity comparison between primer sets
#'
#' @param div1 Diversity dataframe from PS1
#' @param div2 Diversity dataframe from PS2
#' @param group_var Character, grouping variable name
#' @param title Character, plot title
#' @return ggplot object
plot_comparative_shannon <- function(div1, div2, group_var, title) {
  # Add primer set labels
  div1$PrimerSet <- "PS1"
  div2$PrimerSet <- "PS2"

  # Combine
  div_combined <- bind_rows(div1, div2)

  # Plot
  p <- ggplot(div_combined, aes(
    x = .data[[group_var]], y = Shannon,
    fill = PrimerSet
  )) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.7) +
    scale_fill_manual(values = c("PS1" = ps1_color, "PS2" = ps2_color)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = title, x = group_var, y = "Shannon Diversity Index",
      fill = "Primer Set"
    )

  return(p)
}


#' Plot faceted PCoA comparison
#'
#' @param beta1 Beta diversity result from PS1
#' @param beta2 Beta diversity result from PS2
#' @param color_var Character, variable to color points by
#' @param shape_var Character, variable for point shapes (optional)
#' @param title Character, plot title
#' @return ggplot object
plot_faceted_pcoa <- function(beta1, beta2, color_var, shape_var = NULL,
                              title = "PCoA Comparison") {
  # Create dataframes
  pcoa1_df <- data.frame(
    PCoA1 = beta1$pcoa$points[, 1],
    PCoA2 = beta1$pcoa$points[, 2],
    beta1$metadata,
    PrimerSet = "Primerset 1"
  )

  pcoa2_df <- data.frame(
    PCoA1 = beta2$pcoa$points[, 1],
    PCoA2 = beta2$pcoa$points[, 2],
    beta2$metadata,
    PrimerSet = "Primerset 2"
  )

  # Combine
  pcoa_combined <- bind_rows(pcoa1_df, pcoa2_df)

  # Build plot
  if (is.null(shape_var)) {
    p <- ggplot(pcoa_combined, aes(
      x = PCoA1, y = PCoA2,
      color = .data[[color_var]]
    )) +
      geom_point(size = 3, alpha = 0.7)
  } else {
    p <- ggplot(pcoa_combined, aes(
      x = PCoA1, y = PCoA2,
      color = .data[[color_var]],
      shape = .data[[shape_var]]
    )) +
      geom_point(size = 3, alpha = 0.7)
  }

  p <- p +
    facet_wrap(~PrimerSet, scales = "free") +
    theme_minimal() +
    labs(title = title, color = color_var, shape = shape_var)

  return(p)
}


# ===== CROSS-VALIDATION STATISTICS =====

#' Correlate Shannon diversity for shared samples
#'
#' @param div1 Diversity dataframe from PS1
#' @param div2 Diversity dataframe from PS2
#' @return List with correlation test and shared sample count
correlate_shared_samples <- function(div1, div2) {
  # Find shared samples
  shared <- intersect(div1$Sample_Name, div2$Sample_Name)

  if (length(shared) < 3) {
    return(list(
      cor_test = NULL, n_shared = length(shared),
      message = "Insufficient shared samples"
    ))
  }

  # Subset and match order
  div1_shared <- div1[match(shared, div1$Sample_Name), ]
  div2_shared <- div2[match(shared, div2$Sample_Name), ]

  # Pearson correlation
  cor_test <- cor.test(div1_shared$Shannon, div2_shared$Shannon,
    method = "pearson"
  )

  return(list(cor_test = cor_test, n_shared = length(shared)))
}


#' Run Mantel test for beta diversity concordance
#'
#' @param dist1 Distance matrix from PS1
#' @param dist2 Distance matrix from PS2
#' @return Mantel test result
run_mantel_test <- function(dist1, dist2) {
  # Ensure labels match
  if (!identical(labels(dist1), labels(dist2))) {
    warning("Distance matrix labels don't match - attempting to align")
    shared_labels <- intersect(labels(dist1), labels(dist2))

    # Subset to shared samples (requires converting to matrix)
    mat1 <- as.matrix(dist1)
    mat2 <- as.matrix(dist2)
    mat1 <- mat1[shared_labels, shared_labels]
    mat2 <- mat2[shared_labels, shared_labels]
    dist1 <- as.dist(mat1)
    dist2 <- as.dist(mat2)
  }

  # Mantel test
  mantel_result <- mantel(dist1, dist2,
    method = "pearson",
    permutations = 999
  )

  return(mantel_result)
}


#' Check if two p-values agree on significance
#'
#' @param p1 P-value from test 1
#' @param p2 P-value from test 2
#' @param alpha Significance threshold (default 0.05)
#' @return Character: "Yes", "No", or "Partial"
calculate_concordance <- function(p1, p2, alpha = 0.05) {
  if (is.na(p1) || is.na(p2)) {
    return(NA)
  }

  sig1 <- p1 < alpha
  sig2 <- p2 < alpha

  if (sig1 == sig2) {
    return("Yes")
  } else {
    return("No")
  }
}


#' Compare effect sizes between two analyses
#'
#' @param result1 Statistical result from PS1
#' @param result2 Statistical result from PS2
#' @param metric Character: "R2", "eta2", or "F"
#' @return Numeric difference in effect sizes
compare_effect_sizes <- function(result1, result2, metric = "R2") {
  if (metric == "R2") {
    # Extract R² from PERMANOVA
    r2_1 <- result1$R2[1]
    r2_2 <- result2$R2[1]
    return(abs(r2_1 - r2_2))
  } else if (metric == "eta2") {
    # Calculate eta² from ANOVA summary
    ss1 <- summary(result1)[[1]]
    eta2_1 <- ss1$`Sum Sq`[1] / sum(ss1$`Sum Sq`)

    ss2 <- summary(result2)[[1]]
    eta2_2 <- ss2$`Sum Sq`[1] / sum(ss2$`Sum Sq`)

    return(abs(eta2_1 - eta2_2))
  } else if (metric == "F") {
    # Compare F-statistics
    f1 <- result1$F[1]
    f2 <- result2$F[1]
    return(abs(f1 - f2))
  }
}


#' Filter species with low sample counts
#'
#' @param ps_obj phyloseq object
#' @param min_n Minimum number of samples per species (default 3)
#' @return Filtered phyloseq object (or NULL if no species remain)
filter_low_n_species <- function(ps_obj, min_n = 3) {
  meta <- data.frame(sample_data(ps_obj), stringsAsFactors = FALSE)

  # CRITICAL: Use rownames from metadata which match phyloseq sample names

  species_counts <- table(meta$Species)
  species_to_keep <- names(species_counts)[species_counts >= min_n]

  cat("Removing species with < ", min_n, " samples:\n", sep = "")
  removed <- names(species_counts)[species_counts < min_n]
  if (length(removed) > 0) {
    for (sp in removed) {
      cat("  -", sp, "(n =", species_counts[sp], ")\n")
    }
  } else {
    cat("  (none)\n")
  }

  # Check if any species remain
  if (length(species_to_keep) == 0) {
    cat("\n⚠ WARNING: No species have >= ", min_n, " samples!\n", sep = "")
    cat("Cannot proceed with filtering. Returning NULL.\n")
    return(NULL)
  }

  # Get sample IDs to keep (using rownames which match sample_names(ps_obj))
  # Directly filter based on Species column in metadata
  s2k <- rownames(meta)[meta$Species %in% species_to_keep]

  # Safety check
  if (length(s2k) == 0) {
    cat("\n⚠ WARNING: No samples remain after filtering!\n")
    return(NULL)
  }

  # Prune the phyloseq object
  ps_filt <- prune_samples(s2k, ps_obj)

  # CRITICAL FIX: Explicitly rebuild sample_data with correct rownames
  # This ensures metadata rownames match sample_names after all our manipulations
  meta_filt <- data.frame(sample_data(ps_filt), stringsAsFactors = FALSE)
  rownames(meta_filt) <- sample_names(ps_filt)
  sample_data(ps_filt) <- sample_data(meta_filt)

  cat(
    "Retained", length(species_to_keep), "species,",
    nsamples(ps_filt), "samples\n"
  )

  return(ps_filt)
}

# ===== TAXONOMIC HARMONIZATION =====

#' Clean and standardize family-level taxonomy
#'
#' @param family_string Character, raw family name
#' @return Character, cleaned family name
clean_family_name <- function(family_string) {
  if (is.na(family_string) || family_string == "" ||
    grepl("Unknown|unknown|Unassigned|unassigned", family_string)) {
    return(NA_character_)
  }

  # Remove QIIME2 prefix pattern (e.g., "f__", "g__", etc.)
  family_clean <- sub("^[a-z]__", "", family_string)

  # If that didn't work, try more explicit pattern
  family_clean <- sub("^f__", "", family_clean)

  # Remove numerical suffixes (e.g., "_7149")
  family_clean <- sub("_[0-9]+$", "", family_clean)

  # Remove extra whitespace
  family_clean <- trimws(family_clean)

  # Return NA if cleaning resulted in empty string
  if (family_clean == "") {
    return(NA_character_)
  }

  return(family_clean)
}


#' Clean and standardize species names (vectorized)
#'
#' @param species_string Character vector of raw species names
#' @return Character vector of cleaned species names
clean_species_name <- function(species_string) {
  # Convert to character (in case it's a factor)
  result <- as.character(species_string)

  # Handle specific case: P. kingsleyae is Paramormyrops, not Petrocephalus
  result <- gsub("^P\\. kingsleyae", "Paramormyrops kingsleyae", result)

  # Expand other abbreviated genus names
  result <- gsub("^M\\. ", "Marcusenius ", result)
  result <- gsub("^P\\. ", "Paramormyrops ", result)

  # Fix known typos
  result <- gsub("zanchlirostrus", "zanchlirostris", result)

  return(result)
}


#' Extract family-level aggregated counts from phyloseq object
#'
#' @param ps phyloseq object
#' @return Data frame with Sample, Family, Count
extract_family_counts <- function(ps) {
  # Get taxonomy table
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

  # Get OTU table (features as rows, samples as columns)
  otu_df <- as.data.frame(otu_table(ps), stringsAsFactors = FALSE)

  # Ensure they align
  if (taxa_are_rows(ps)) {
    # Already correct orientation
  } else {
    otu_df <- t(otu_df)
  }

  # Add cleaned family names to OTU table
  otu_df$Family_clean <- sapply(tax_df[rownames(otu_df), "Family"], clean_family_name)

  # Remove ASVs without family assignment
  otu_df <- otu_df %>% filter(!is.na(Family_clean))

  # Reshape to long format
  otu_long <- otu_df %>%
    pivot_longer(cols = -Family_clean, names_to = "SampleID", values_to = "Count") %>%
    filter(Count > 0)

  # Aggregate to family level per sample
  family_counts <- otu_long %>%
    group_by(SampleID, Family_clean) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    rename(Family = Family_clean)

  return(family_counts)
}


#' Create global color palette for family-level taxa
#'
#' @param ps1 phyloseq object for primer set 1
#' @param ps2 phyloseq object for primer set 2
#' @param top_n Integer, number of top families to include (default 25)
#' @return Named character vector mapping family names to colors
create_global_family_palette <- function(ps1, ps2, min_families = 50) {
  # Transform to relative abundance
  ps1_rel <- transform_sample_counts(ps1, function(x) x / sum(x))
  ps2_rel <- transform_sample_counts(ps2, function(x) x / sum(x))

  # Get family-level abundance from both datasets
  ps1_fam_agg <- tax_glom(ps1_rel, taxrank = "Family", NArm = FALSE)
  ps2_fam_agg <- tax_glom(ps2_rel, taxrank = "Family", NArm = FALSE)

  # Extract mean abundances
  ps1_abund <- taxa_sums(ps1_fam_agg) / nsamples(ps1_fam_agg)
  ps2_abund <- taxa_sums(ps2_fam_agg) / nsamples(ps2_fam_agg)

  # Get family names
  ps1_tax <- tax_table(ps1_fam_agg)[, "Family"] %>%
    as.vector() %>%
    sapply(clean_family_name)
  ps2_tax <- tax_table(ps2_fam_agg)[, "Family"] %>%
    as.vector() %>%
    sapply(clean_family_name)

  # Combine and aggregate
  all_abund <- data.frame(
    Family = c(ps1_tax, ps2_tax),
    Abundance = c(ps1_abund, ps2_abund),
    stringsAsFactors = FALSE
  ) %>%
    group_by(Family) %>%
    summarise(MeanAbund = mean(Abundance), .groups = "drop") %>%
    filter(!is.na(Family) & Family != "") %>%
    arrange(desc(MeanAbund))

  n_families <- nrow(all_abund)

  # Generate highly distinct colors using multiple strategies
  # 1. Start with hand-picked distinct colors for most common families
  # 2. Supplement with algorithmically spaced colors for remaining families

  # Set of maximally distinct base colors (avoid greys - reserved for "Other")
  # Expanded palette with more distinct hues
  distinct_colors <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377",
    "#EE3377", "#009988", "#CC3311", "#EE7733", "#0077BB", "#33BBEE",
    "#CC6677", "#882255", "#AA4499", "#117733", "#332288", "#44AA99", "#999933",
    "#DDCC77", "#661100", "#6699CC", "#AA4466", "#771122", "#DD7788",
    "#88CCEE", "#332288", "#117733", "#999933", "#882255", "#CC6677"
  )

  if (n_families <= length(distinct_colors)) {
    # Use distinct colors directly
    family_colors <- setNames(
      distinct_colors[1:n_families],
      all_abund$Family
    )
  } else {
    # For many families, use distinct colors for top families
    # then fill in with algorithmically spaced colors
    n_distinct <- length(distinct_colors)
    n_additional <- n_families - n_distinct

    # Generate additional colors using rainbow with maximum hue separation
    additional_colors <- grDevices::rainbow(n_additional, s = 0.8, v = 0.8)

    # Combine: distinct colors first (for most abundant), then additional
    all_colors <- c(distinct_colors, additional_colors)

    family_colors <- setNames(
      all_colors[1:n_families],
      all_abund$Family
    )
  }

  # Add "Other" category in gray
  family_colors["Other"] <- "#808080"

  return(family_colors)
}


# ===== SELECTIVITY CALCULATIONS =====

#' Calculate Ivlev's Electivity Index
#'
#' @param r Numeric, proportion in diet (use)
#' @param p Numeric, proportion in environment (availability)
#' @return Numeric, Ivlev's E (-1 to +1)
ivlev_index <- function(r, p) {
  if (r + p == 0) {
    return(NA_real_)
  }
  return((r - p) / (r + p))
}


#' Calculate Chesson's Alpha
#'
#' @param r Numeric vector, proportions in diet
#' @param p Numeric vector, proportions in environment
#' @return Numeric vector, Chesson's alpha values
chesson_alpha <- function(r, p) {
  # Avoid division by zero
  ratio <- ifelse(p > 0, r / p, NA_real_)
  sum_ratio <- sum(ratio, na.rm = TRUE)

  if (sum_ratio == 0) {
    return(rep(NA_real_, length(r)))
  }

  alpha <- ratio / sum_ratio
  return(alpha)
}


#' Calculate Jacobs' Index
#'
#' @param r Numeric, proportion in diet
#' @param p Numeric, proportion in environment
#' @return Numeric, Jacobs' D (-1 to +1)
jacobs_index <- function(r, p) {
  denominator <- r + p - 2 * r * p
  if (denominator == 0) {
    return(NA_real_)
  }
  return((r - p) / denominator)
}


#' Calculate all selectivity indices for a fish group
#'
#' @param diet_counts Data frame with Family and Count columns (diet)
#' @param env_counts Data frame with Family and Count columns (environment)
#' @return Data frame with Family, Availability_prop, Use_prop, and selectivity indices
calculate_selectivity <- function(diet_counts, env_counts) {
  # Calculate proportions
  diet_prop <- diet_counts %>%
    mutate(Use_prop = Count / sum(Count)) %>%
    select(Family, Use_prop)

  env_prop <- env_counts %>%
    mutate(Availability_prop = Count / sum(Count)) %>%
    select(Family, Availability_prop)

  # Get all families present in either dataset
  all_families <- unique(c(diet_prop$Family, env_prop$Family))

  # Create complete dataframe
  selectivity_df <- data.frame(Family = all_families) %>%
    left_join(env_prop, by = "Family") %>%
    left_join(diet_prop, by = "Family") %>%
    replace_na(list(Availability_prop = 0, Use_prop = 0))

  # Calculate Chesson's alpha for all families
  chesson_vals <- chesson_alpha(selectivity_df$Use_prop, selectivity_df$Availability_prop)

  # Calculate indices
  selectivity_df <- selectivity_df %>%
    mutate(
      Ivlev_E = mapply(ivlev_index, Use_prop, Availability_prop),
      Jacobs_D = mapply(jacobs_index, Use_prop, Availability_prop),
      Chesson_Alpha = chesson_vals,
      Selection_Type = case_when(
        is.na(Ivlev_E) ~ "Undefined",
        Ivlev_E > 0.2 ~ "Preference",
        Ivlev_E < -0.2 ~ "Avoidance",
        TRUE ~ "Neutral"
      )
    )

  return(selectivity_df)
}


# ===== VISUALIZATION =====

#' Create heatmap of SIMPER analysis results
#'
#' @param simper_result SIMPER result object from vegan::simper()
#' @param ps phyloseq object
#' @param top_n Integer, number of top families to show per comparison (default 10)
#' @param primerset_name Character, name for plot title
#' @return ggplot object
plot_simper_heatmap <- function(simper_result, ps, top_n = 10,
                                primerset_name = "") {
  # Extract pairwise comparisons
  comparison_list <- list()

  for (comp_name in names(simper_result)) {
    comp_summary <- summary(simper_result)[[comp_name]]

    if (!is.data.frame(comp_summary)) {
      comp_summary <- as.data.frame(comp_summary)
    }

    if (nrow(comp_summary) > 0) {
      top_asvs <- head(rownames(comp_summary), top_n)

      # Get family names from taxonomy
      tax_df <- as.data.frame(tax_table(ps))
      families <- sapply(top_asvs, function(asv) {
        if (asv %in% rownames(tax_df)) {
          fam <- tax_df[asv, "Family"]
          if (is.na(fam) || fam == "") {
            return(tax_df[asv, "Order"])
          }
          return(clean_family_name(fam))
        }
        return("Unknown")
      })

      # Aggregate by family
      for (i in seq_along(families)) {
        comparison_list[[length(comparison_list) + 1]] <- data.frame(
          Comparison = gsub("_", " vs ", comp_name),
          Family = families[i],
          Contribution = comp_summary[i, "average"],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  all_data <- do.call(rbind, comparison_list)

  # Get top families overall
  family_totals <- all_data %>%
    group_by(Family) %>%
    summarise(Total = sum(Contribution), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = top_n)

  plot_data <- all_data %>% filter(Family %in% family_totals$Family)

  # Create heatmap
  p <- ggplot(plot_data, aes(x = Comparison, y = Family, fill = Contribution)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_viridis_c(
      option = "plasma",
      name = "Contribution\nto Dissimilarity"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0(
        "Top Prey Contributing to Species Differences (",
        primerset_name, ")"
      ),
      x = "Species Comparison", y = "Prey Family"
    )

  return(p)
}


#' Create heatmap of selectivity indices
#'
#' @param selectivity_data Data frame with columns: Fish_Group, Family, Ivlev_E
#' @param title Character, plot title
#' @return ggplot object
plot_selectivity_heatmap <- function(selectivity_data, title = "Prey Selectivity") {
  # Filter to families present in at least 2 fish groups
  family_counts <- selectivity_data %>%
    filter(!is.na(Ivlev_E)) %>%
    group_by(Family) %>%
    summarise(n_groups = n_distinct(Fish_Group)) %>%
    filter(n_groups >= 2)

  plot_data <- selectivity_data %>%
    filter(Family %in% family_counts$Family)

  p <- ggplot(plot_data, aes(x = Family, y = Fish_Group, fill = Ivlev_E)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(
      low = avoidance_color, mid = neutral_color,
      high = preference_color, midpoint = 0,
      limits = c(-1, 1), na.value = "gray90",
      name = "Ivlev's E"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(title = title, x = "Prey Family", y = "Fish Group")

  return(p)
}


#' Create barplot comparing availability vs use
#'
#' @param availability_data Data frame with Site, Family, Proportion
#' @param use_data Data frame with Site, Fish_Group, Family, Proportion
#' @param site_filter Character vector, sites to include
#' @return ggplot object
plot_availability_vs_use <- function(availability_data, use_data, site_filter = NULL) {
  if (!is.null(site_filter)) {
    availability_data <- availability_data %>% filter(Site %in% site_filter)
    use_data <- use_data %>% filter(Site %in% site_filter)
  }

  # Prepare availability data
  avail_plot <- availability_data %>%
    mutate(Type = "Available")

  # Prepare use data (aggregate across fish groups per site)
  use_plot <- use_data %>%
    group_by(Site, Family) %>%
    summarise(Proportion = mean(Proportion), .groups = "drop") %>%
    mutate(Type = "Used")

  # Combine
  combined <- bind_rows(avail_plot, use_plot)

  # Plot
  p <- ggplot(combined, aes(x = Type, y = Proportion, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Site, ncol = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    labs(
      title = "Prey Availability vs. Use by Site",
      x = "", y = "Proportion", fill = "Prey Family"
    ) +
    scale_fill_viridis_d()

  return(p)
}


# ===== FUZZY MATCHING FOR PHYLOGENY =====

#' Parse tree tip labels to extract genus and species
#'
#' @param tip_labels Character vector of tree tip labels
#' @return Character vector of "Genus species" format
parse_tree_tips <- function(tip_labels) {
  # Split on underscores and extract first two parts
  parsed <- sapply(tip_labels, function(tip) {
    parts <- strsplit(tip, "_")[[1]]
    if (length(parts) >= 2) {
      paste(parts[1], parts[2])
    } else {
      tip # Return original if can't parse
    }
  }, USE.NAMES = FALSE)

  return(parsed)
}


#' Normalize Bale Creek species names using genus context
#'
#' @param species_vec Character vector of species names (may be abbreviated)
#' @param genus_vec Character vector of genus names from metadata
#' @return Character vector of normalized "Genus species" format
normalize_bale_species <- function(species_vec, genus_vec) {
  normalized <- character(length(species_vec))

  for (i in seq_along(species_vec)) {
    sp <- species_vec[i]
    genus <- genus_vec[i]

    # Check if species is abbreviated (starts with single letter + period)
    if (grepl("^[A-Z]\\. ", sp)) {
      # Replace abbreviated genus with full genus from metadata
      normalized[i] <- sub("^[A-Z]\\. ", paste0(genus, " "), sp)
    } else {
      # Already in full format
      normalized[i] <- sp
    }
  }

  return(normalized)
}


#' Fuzzy match species names with Levenshtein distance
#'
#' @param query Character, species name to match
#' @param targets Character vector, candidate species names
#' @param max_dist Integer, maximum Levenshtein distance (default 2)
#' @return List with matched_target, distance, match_type
fuzzy_match_species <- function(query, targets, max_dist = 2) {
  # Try exact match first
  exact_match <- which(targets == query)
  if (length(exact_match) > 0) {
    return(list(
      matched_target = targets[exact_match[1]],
      distance = 0,
      match_type = "exact"
    ))
  }

  # Calculate Levenshtein distances
  distances <- stringdist::stringdist(query, targets, method = "lv")

  # Find best match within threshold
  best_idx <- which.min(distances)
  best_dist <- distances[best_idx]

  if (best_dist <= max_dist) {
    return(list(
      matched_target = targets[best_idx],
      distance = best_dist,
      match_type = "fuzzy"
    ))
  }

  # No match found
  return(list(
    matched_target = NA_character_,
    distance = NA_integer_,
    match_type = "no_match"
  ))
}


#' Match Bale Creek species to phylogeny with fuzzy matching
#'
#' @param ps_bale phyloseq object for Bale Creek samples
#' @param phylo_tree phylo object (tree)
#' @return List with matched_species, match_table, n_matched, n_species
match_bale_species_to_phylogeny <- function(ps_bale, phylo_tree) {
  # Define special case dictionary
  # Map full species names (as they appear in metadata) to tree tip PATTERNS
  # Note: tree tips have suffixes, so we use patterns like "^Paramormyrops_SZA"
  special_cases <- list(
    "Paramormyrops spp. 'SZA'" = "^Paramormyrops_SZA",
    "Paramormyrops spp. SZA" = "^Paramormyrops_SZA",
    "Paramormyrops spp. 'MAG'" = "^Paramormyrops_mag",
    "Paramormyrops spp. MAG" = "^Paramormyrops_mag",
    "Paramormyrops spp. 'TEN'" = "^Paramormyrops_TEN",
    "Paramormyrops spp. TEN" = "^Paramormyrops_TEN",
    "Paramormyrops spp. 'TEU'" = "^Paramormyrops_TEU",
    "Paramormyrops spp. TEU" = "^Paramormyrops_TEU",
    "Paramormyrops spp. 'VAD'" = "^Paramormyrops_batesii",
    "Paramormyrops spp. VAD" = "^Paramormyrops_batesii",
    "Marcusenius moori" = "^Marcusenius_moorii",
    "Mormyrops zanchlirostris" = "^Mormyrops_zanclirostris",
    "Petrocephalus simus" = "^Petrocephalus_balayi"
  )

  # Get Bale species from phyloseq metadata
  meta <- data.frame(sample_data(ps_bale), stringsAsFactors = FALSE)

  # Ensure Species and Genus columns are character, not factor
  meta$Species <- as.character(meta$Species)
  meta$Genus <- as.character(meta$Genus)

  # Get unique species
  bale_species_raw <- unique(meta$Species)

  # Create genus lookup for each species
  genus_lookup <- character(length(bale_species_raw))
  names(genus_lookup) <- bale_species_raw
  for (sp in bale_species_raw) {
    genus_lookup[sp] <- meta$Genus[meta$Species == sp][1]
  }

  # Normalize Bale species names (expand abbreviations if needed)
  bale_species_normalized <- character(length(bale_species_raw))
  for (i in seq_along(bale_species_raw)) {
    sp <- bale_species_raw[i]
    genus <- genus_lookup[sp]

    # Check if species is abbreviated (starts with single letter + period)
    if (grepl("^[A-Z]\\. ", sp)) {
      # Replace abbreviated genus with full genus from metadata
      bale_species_normalized[i] <- sub("^[A-Z]\\. ", paste0(genus, " "), sp)
    } else {
      # Already in full format
      bale_species_normalized[i] <- sp
    }
  }

  # Parse tree tips
  tree_tips <- phylo_tree$tip.label
  tree_species <- parse_tree_tips(tree_tips)
  names(tree_species) <- tree_tips

  # Debug output
  cat("\nDEBUG - Bale species to match:\n")
  for (i in seq_along(bale_species_raw)) {
    cat("  ", i, ": ", bale_species_raw[i], " --> ", bale_species_normalized[i], "\n", sep = "")
  }

  # Prepare results table
  match_table <- data.frame(
    bale_species = bale_species_raw,
    normalized_species = bale_species_normalized,
    tree_tip = NA_character_,
    tree_species = NA_character_,
    match_type = NA_character_,
    distance = NA_integer_,
    stringsAsFactors = FALSE
  )

  # Match each Bale species
  for (i in seq_along(bale_species_normalized)) {
    query <- bale_species_normalized[i]

    # Check special cases first (try both with and without quotes)
    query_clean <- gsub("'", "", query) # Remove quotes for matching

    # Try to find matching special case
    special_pattern <- NULL
    if (query %in% names(special_cases)) {
      special_pattern <- special_cases[[query]]
    } else if (query_clean %in% names(special_cases)) {
      special_pattern <- special_cases[[query_clean]]
    }

    if (!is.null(special_pattern)) {
      # Find tree tip matching this pattern
      matched_tips <- grep(special_pattern, tree_tips, value = TRUE)

      if (length(matched_tips) > 0) {
        # Take first match
        match_table$tree_tip[i] <- matched_tips[1]
        match_table$tree_species[i] <- tree_species[matched_tips[1]]
        match_table$match_type[i] <- "special_case"
        match_table$distance[i] <- NA_integer_
        next
      }
    }

    # Try fuzzy matching on parsed species names
    match_result <- fuzzy_match_species(query, tree_species, max_dist = 2)

    if (!is.na(match_result$matched_target)) {
      # Find corresponding tip label
      matched_tip <- names(tree_species)[tree_species == match_result$matched_target][1]

      match_table$tree_tip[i] <- matched_tip
      match_table$tree_species[i] <- match_result$matched_target
      match_table$match_type[i] <- match_result$match_type
      match_table$distance[i] <- match_result$distance
    }
  }

  # Count successful matches
  n_matched <- sum(!is.na(match_table$tree_tip))
  n_species <- nrow(match_table)

  # Extract matched tips for downstream use
  matched_tips <- match_table$tree_tip[!is.na(match_table$tree_tip)]

  return(list(
    matched_species = matched_tips,
    match_table = match_table,
    n_matched = n_matched,
    n_species = n_species
  ))
}


# ===== ADDITIONAL ANALYSIS FUNCTIONS =====

#' Aggregate phyloseq object to family level
#'
#' @param ps phyloseq object
#' @return phyloseq object aggregated to family level
aggregate_to_family <- function(ps) {
  ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)
  return(ps_family)
}


#' Create stacked bar chart of prey composition at family level
#'
#' @param ps phyloseq object (should be aggregated to family level)
#' @param group_var Character, grouping variable for x-axis
#' @param top_n Integer, number of top families to show (rest as "Other")
#' @param title Character, plot title
#' @return ggplot object
plot_composition_bars <- function(ps, group_var, top_n = 15, title = "") {
  # Melt to long format (before normalization)
  ps_melt <- psmelt(ps)

  # Clean family names
  ps_melt$Family_clean <- sapply(ps_melt$Family, clean_family_name)

  # Aggregate by group_var and Family, summing counts
  ps_agg <- ps_melt %>%
    group_by(.data[[group_var]], Family_clean) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # Normalize within each group to get relative abundance
  ps_agg <- ps_agg %>%
    group_by(.data[[group_var]]) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>%
    ungroup()

  # Identify top N families by mean relative abundance across groups
  top_families <- ps_agg %>%
    group_by(Family_clean) %>%
    summarise(mean_abund = mean(RelAbund), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(Family_clean)

  # Categorize families and aggregate "Other"
  ps_final <- ps_agg %>%
    mutate(Family_display = ifelse(Family_clean %in% top_families,
      Family_clean, "Other"
    )) %>%
    group_by(.data[[group_var]], Family_display) %>%
    summarise(RelAbund = sum(RelAbund), .groups = "drop")

  # Reorder for plotting (top families first, then Other)
  ps_final$Family_display <- factor(ps_final$Family_display,
    levels = c(top_families, "Other")
  )

  # Create plot with consistent color scale
  p <- ggplot(ps_final, aes(
    x = .data[[group_var]], y = RelAbund,
    fill = Family_display
  )) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = GLOBAL_FAMILY_COLORS,
      name = "Prey Family",
      na.value = "#808080",
      drop = FALSE
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.key.size = unit(0.35, "cm"),
      legend.box = "vertical"
    ) +
    guides(fill = guide_legend(
      ncol = 4, title.position = "top",
      byrow = TRUE, nrow = 3
    )) +
    labs(title = title, x = group_var, y = "Relative Abundance") +
    scale_y_continuous(labels = scales::percent)

  return(p)
}


#' Variance partitioning for PERMANOVA
#'
#' @param beta_result Beta diversity result from run_beta_diversity_analysis()
#' @param var1_formula Character, formula for first variable (e.g., "~ Region")
#' @param var2_formula Character, formula for second variable (e.g., "~ Waveform_Type")
#' @return List with varpart result and sequential PERMANOVA results
variance_partition_permanova <- function(beta_result, var1_formula, var2_formula) {
  dist_mat <- beta_result$distance
  meta <- beta_result$metadata

  # Method 1: varpart
  varpart_result <- varpart(dist_mat,
    as.formula(var1_formula),
    as.formula(var2_formula),
    data = meta
  )

  # Method 2: Sequential PERMANOVA (variable 1 first)
  formula_seq1 <- paste(
    "dist_mat", var1_formula, "+",
    gsub("~\\s*", "", var2_formula)
  )
  perm_seq1 <- adonis2(as.formula(formula_seq1),
    data = meta, permutations = 999, by = "terms"
  )

  # Method 3: Sequential PERMANOVA (variable 2 first)
  formula_seq2 <- paste(
    "dist_mat", var2_formula, "+",
    gsub("~\\s*", "", var1_formula)
  )
  perm_seq2 <- adonis2(as.formula(formula_seq2),
    data = meta, permutations = 999, by = "terms"
  )

  # Method 4: Marginal effects
  formula_marginal <- paste(
    "dist_mat", var1_formula, "+",
    gsub("~\\s*", "", var2_formula)
  )
  perm_marginal <- adonis2(as.formula(formula_marginal),
    data = meta, permutations = 999, by = "margin"
  )

  return(list(
    varpart = varpart_result,
    sequential_var1_first = perm_seq1,
    sequential_var2_first = perm_seq2,
    marginal = perm_marginal
  ))
}


#' Calculate diet overlap using Schoener and Pianka indices
#'
#' @param ps phyloseq object
#' @param group_var Character, grouping variable (e.g., "Species")
#' @return List with Schoener and Pianka matrices, species_list, and mean_abundances
calculate_diet_overlap <- function(ps, group_var) {
  # Convert to relative abundance
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

  # Get OTU table (taxa as rows, samples as columns)
  otu_mat <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) {
    otu_mat <- t(otu_mat)
  }

  # Get metadata
  meta <- data.frame(sample_data(ps_rel))

  # Calculate mean abundance per group
  groups <- unique(meta[[group_var]])
  mean_abundances <- matrix(0, nrow = ntaxa(ps_rel), ncol = length(groups))
  rownames(mean_abundances) <- taxa_names(ps_rel)
  colnames(mean_abundances) <- groups

  for (g in groups) {
    samples_in_group <- rownames(meta)[meta[[group_var]] == g]
    if (length(samples_in_group) == 1) {
      mean_abundances[, g] <- otu_mat[, samples_in_group]
    } else {
      mean_abundances[, g] <- rowMeans(otu_mat[, samples_in_group, drop = FALSE])
    }
  }

  # Capture species names from the mean_abundances matrix
  # mean_abundances has: rows = taxa, columns = groups (species)
  species_names <- colnames(mean_abundances)

  # Calculate overlap indices
  # NOTE: niche.overlap calculates overlap between COLUMNS
  # mean_abundances has groups/species in columns, which is what we want
  schoener <- niche.overlap(mean_abundances, method = "schoener")
  pianka <- niche.overlap(mean_abundances, method = "pianka")

  # CRITICAL: niche.overlap returns a dist object, convert to matrix
  # dist objects have diagonal = 0 (for distances), but overlap indices should have diagonal = 1
  schoener <- as.matrix(schoener)
  pianka <- as.matrix(pianka)

  # Set diagonal to 1 (perfect overlap with self)
  diag(schoener) <- 1
  diag(pianka) <- 1

  # Assign species names as dimnames
  rownames(schoener) <- species_names
  colnames(schoener) <- species_names
  rownames(pianka) <- species_names
  colnames(pianka) <- species_names

  return(list(
    schoener = schoener,
    pianka = pianka,
    mean_abundances = mean_abundances,
    species_list = species_names
  ))
}


# ===== TURNOVER vs ABUNDANCE ANALYSIS =====

#' Calculate ASV overlap between two groups using Jaccard similarity
#'
#' @param ps phyloseq object
#' @param group_var Character, metadata column name (e.g., "Sex")
#' @param group1 Character, first group name (e.g., "Male")
#' @param group2 Character, second group name (e.g., "Female")
#' @return List with overlap metrics and ASV vectors
calculate_asv_overlap_jaccard <- function(ps, group_var, group1, group2) {
  # Get metadata
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

  # Get OTU table (taxa as rows)
  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get samples for each group
  samples_g1 <- rownames(meta)[meta[[group_var]] == group1]
  samples_g2 <- rownames(meta)[meta[[group_var]] == group2]

  # Get ASVs present in each group (>0 reads in at least one sample)
  asvs_g1 <- rownames(otu_mat)[rowSums(otu_mat[, samples_g1, drop = FALSE]) > 0]
  asvs_g2 <- rownames(otu_mat)[rowSums(otu_mat[, samples_g2, drop = FALSE]) > 0]

  # Set operations
  shared_asvs <- intersect(asvs_g1, asvs_g2)
  unique_g1 <- setdiff(asvs_g1, asvs_g2)
  unique_g2 <- setdiff(asvs_g2, asvs_g1)

  # Jaccard similarity
  union_asvs <- union(asvs_g1, asvs_g2)
  jaccard <- length(shared_asvs) / length(union_asvs)

  return(list(
    n_shared = length(shared_asvs),
    n_unique_g1 = length(unique_g1),
    n_unique_g2 = length(unique_g2),
    n_total = length(union_asvs),
    jaccard = jaccard,
    shared_asvs = shared_asvs,
    unique_g1 = unique_g1,
    unique_g2 = unique_g2
  ))
}


#' Test abundance differences for shared ASVs between two groups
#'
#' @param ps phyloseq object
#' @param shared_asvs Character vector of ASV IDs to test
#' @param group_var Character, grouping variable (e.g., "Sex")
#' @param group1 Character, first group name
#' @param group2 Character, second group name
#' @return Data frame with test results for each ASV
test_shared_asv_abundances <- function(ps, shared_asvs, group_var, group1, group2) {
  # Get metadata
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

  # Get OTU table (taxa as rows)
  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get taxonomy
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

  # Get samples for each group (ensure they match OTU matrix column names)
  samples_g1 <- rownames(meta)[meta[[group_var]] == group1]
  samples_g2 <- rownames(meta)[meta[[group_var]] == group2]

  # Ensure samples exist in OTU matrix
  samples_g1 <- intersect(samples_g1, colnames(otu_mat))
  samples_g2 <- intersect(samples_g2, colnames(otu_mat))

  # Initialize results dataframe
  results <- data.frame(
    ASV = shared_asvs,
    Mean_G1 = NA_real_,
    Mean_G2 = NA_real_,
    Median_G1 = NA_real_,
    Median_G2 = NA_real_,
    Wilcox_W = NA_real_,
    Wilcox_p = NA_real_,
    Family = NA_character_,
    stringsAsFactors = FALSE
  )

  # Test each ASV
  for (i in seq_along(shared_asvs)) {
    asv <- shared_asvs[i]

    # Extract abundances (ensure ASV exists in matrix)
    if (asv %in% rownames(otu_mat)) {
      # Extract as numeric vectors (not otu_table objects)
      g1_abundances <- as.numeric(otu_mat[asv, samples_g1])
      g2_abundances <- as.numeric(otu_mat[asv, samples_g2])

      # Calculate summary stats
      results$Mean_G1[i] <- mean(g1_abundances)
      results$Mean_G2[i] <- mean(g2_abundances)
      results$Median_G1[i] <- median(g1_abundances)
      results$Median_G2[i] <- median(g2_abundances)

      # Wilcoxon test
      wt <- wilcox.test(g1_abundances, g2_abundances, exact = FALSE)
      results$Wilcox_W[i] <- wt$statistic
      results$Wilcox_p[i] <- wt$p.value
    }

    # Get family annotation
    if (asv %in% rownames(tax_df)) {
      results$Family[i] <- tax_df[asv, "Family"]
    }
  }

  # FDR correction
  results$Wilcox_p_adj <- p.adjust(results$Wilcox_p, method = "fdr")
  results$Significant <- results$Wilcox_p_adj < 0.05

  # Sort by p-value
  results <- results[order(results$Wilcox_p), ]

  return(results)
}


#' Partition Bray-Curtis dissimilarity into turnover and abundance components
#'
#' @param ps phyloseq object
#' @param group_var Character, grouping variable
#' @param group1 Character, first group name
#' @param group2 Character, second group name
#' @return List with partitioning metrics
partition_bray_curtis_turnover_abundance <- function(ps, group_var, group1, group2) {
  # Get metadata
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

  # Get OTU table (samples as rows for vegdist)
  otu_mat <- as.matrix(otu_table(ps))
  if (taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get samples for each group
  samples_g1 <- rownames(meta)[meta[[group_var]] == group1]
  samples_g2 <- rownames(meta)[meta[[group_var]] == group2]

  # Calculate BC on full abundance data
  dist_total <- vegdist(otu_mat, method = "bray")

  # Calculate turnover component using Jaccard (designed for presence/absence)
  # Jaccard is more appropriate than BC on binary data
  otu_binary <- (otu_mat > 0) * 1
  dist_turnover <- vegdist(otu_binary, method = "jaccard")

  # Extract between-group distances
  dist_mat_total <- as.matrix(dist_total)
  dist_mat_turnover <- as.matrix(dist_turnover)

  # Get all pairwise distances between groups
  between_total <- dist_mat_total[samples_g1, samples_g2]
  between_turnover <- dist_mat_turnover[samples_g1, samples_g2]

  # Calculate mean dissimilarities
  mean_bc_total <- mean(between_total)
  mean_jaccard_turnover <- mean(between_turnover)

  # For BC partitioning, scale Jaccard to be comparable to BC
  # Use the min of Jaccard and BC_total to ensure turnover <= total
  mean_bc_turnover <- min(mean_jaccard_turnover, mean_bc_total)
  mean_bc_abundance <- mean_bc_total - mean_bc_turnover

  # Calculate proportions
  prop_turnover <- mean_bc_turnover / mean_bc_total
  prop_abundance <- mean_bc_abundance / mean_bc_total

  return(list(
    bc_total = mean_bc_total,
    bc_turnover = mean_bc_turnover,
    bc_abundance = mean_bc_abundance,
    prop_turnover = prop_turnover,
    prop_abundance = prop_abundance,
    n_comparisons = length(between_total)
  ))
}


#' Create heatmap showing shared ASV abundances across samples
#'
#' @param ps phyloseq object
#' @param shared_asvs Character vector of ASVs to plot
#' @param group_var Character, grouping variable (e.g., "Sex")
#' @param test_results Data frame from test_shared_asv_abundances (optional)
#' @param top_n Integer, max ASVs to show (default 30)
#' @param title Character, plot title
#' @return pheatmap object
plot_asv_abundance_heatmap <- function(ps, shared_asvs, group_var,
                                       test_results = NULL, top_n = 30,
                                       title = "Shared ASV Abundances") {
  # Get OTU table (taxa as rows)
  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get taxonomy
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

  # Prioritize ASVs to plot
  if (!is.null(test_results)) {
    # First get significant ASVs
    sig_asvs <- test_results$ASV[test_results$Significant]

    # Then get top abundant non-significant to reach top_n
    mean_abund <- rowMeans(otu_mat[shared_asvs, , drop = FALSE])
    names(mean_abund) <- shared_asvs

    # Remove significant ones from abundance ranking
    nonsig_asvs <- setdiff(shared_asvs, sig_asvs)
    nonsig_abund <- mean_abund[nonsig_asvs]
    nonsig_abund <- sort(nonsig_abund, decreasing = TRUE)

    # Combine: all significant + top abundant non-significant
    n_to_add <- max(0, top_n - length(sig_asvs))
    asvs_to_plot <- c(sig_asvs, names(head(nonsig_abund, n_to_add)))
  } else {
    # Just use top by mean abundance
    mean_abund <- rowMeans(otu_mat[shared_asvs, , drop = FALSE])
    asvs_to_plot <- names(sort(mean_abund, decreasing = TRUE)[1:min(top_n, length(shared_asvs))])
  }

  # Transform to relative abundance
  plot_matrix <- otu_mat[asvs_to_plot, , drop = FALSE]
  plot_matrix <- t(apply(plot_matrix, 1, function(x) x / sum(otu_mat) * 100))

  # Add family names to row labels
  row_labels <- sapply(asvs_to_plot, function(asv) {
    fam <- tax_df[asv, "Family"]
    if (is.na(fam) || fam == "") fam <- "Unknown"
    fam_clean <- clean_family_name(fam)
    paste0(substr(asv, 1, 8), " (", fam_clean, ")")
  })
  rownames(plot_matrix) <- row_labels

  # Create column annotation
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  annotation_col <- data.frame(
    row.names = rownames(meta),
    Group = meta[[group_var]]
  )
  colnames(annotation_col) <- group_var

  # Create row annotation if test results provided
  annotation_row <- NULL
  if (!is.null(test_results)) {
    sig_status <- ifelse(asvs_to_plot %in% test_results$ASV[test_results$Significant],
      "Yes", "No"
    )
    annotation_row <- data.frame(
      row.names = row_labels,
      Significant = sig_status
    )
  }

  # Create heatmap
  pheatmap(plot_matrix,
    scale = "row",
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    main = title,
    fontsize_row = 6,
    fontsize_col = 6,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
}


# ===== SIMPER DRILL-DOWN ANALYSIS =====

#' Subset phyloseq object to Chironomidae (legacy function)
#'
#' @param ps phyloseq object
#' @param primerset_name Character, primer set label for output
#' @return List with ps (subset phyloseq), metadata, n_asvs, or NULL if insufficient data
subset_to_chironomidae <- function(ps, primerset_name = "") {
  cat("\n", primerset_name, " - Subsetting to Chironomidae\n", sep = "")
  result <- subset_to_family(ps, "Chironomidae", min_asvs = 2, min_samples = 5)
  return(result)
}


#' Subset phyloseq object to a specific family
#'
#' @param ps phyloseq object
#' @param family_name Character, family name to subset to (e.g., "Chironomidae")
#' @param min_asvs Integer, minimum ASVs required (default 2)
#' @param min_samples Integer, minimum samples with family reads (default 5)
#' @return List with ps (subset phyloseq), metadata, n_asvs, or NULL if insufficient data
subset_to_family <- function(ps, family_name, min_asvs = 2, min_samples = 5) {
  # Get taxonomy table
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

  # Clean family names in taxonomy
  tax_df$Family_clean <- sapply(tax_df$Family, clean_family_name)

  # Find ASVs in the target family (handle QIIME2 suffixes)
  family_asvs <- rownames(tax_df)[
    !is.na(tax_df$Family_clean) &
      grepl(paste0("^", family_name), tax_df$Family_clean)
  ]

  if (length(family_asvs) < min_asvs) {
    cat(
      "  WARNING: Only", length(family_asvs), "ASVs in", family_name,
      "- need at least", min_asvs, "\n"
    )
    return(NULL)
  }

  # Subset to family ASVs
  ps_family <- prune_taxa(family_asvs, ps)

  # Remove samples with zero reads
  ps_family <- prune_samples(sample_sums(ps_family) > 0, ps_family)

  if (nsamples(ps_family) < min_samples) {
    cat(
      "  WARNING: Only", nsamples(ps_family), "samples with", family_name,
      "- need at least", min_samples, "\n"
    )
    return(NULL)
  }

  # Get metadata
  meta <- data.frame(sample_data(ps_family), stringsAsFactors = FALSE)

  return(list(
    ps = ps_family,
    metadata = meta,
    n_asvs = ntaxa(ps_family)
  ))
}


#' Run complete drill-down analysis for a single family
#'
#' @param ps phyloseq object
#' @param family_name Character, family name to analyze
#' @param group_var Character, grouping variable (e.g., "Sex")
#' @param group1 Character, first group name (e.g., "Male")
#' @param group2 Character, second group name (e.g., "Female")
#' @param primerset_name Character, primer set label for output
#' @return List with all drill-down results, or NULL if insufficient data
run_family_drilldown <- function(ps, family_name, group_var, group1, group2,
                                 primerset_name) {
  cat("\n--- Family:", family_name, "---\n")

  # Subset to family
  family_subset <- subset_to_family(ps, family_name)

  if (is.null(family_subset)) {
    cat("  Skipping", family_name, "- insufficient data\n")
    return(NULL)
  }

  ps_fam <- family_subset$ps

  cat("  ASVs:", family_subset$n_asvs, "\n")
  cat("  Samples:", nsamples(ps_fam), "\n")

  # Check group sizes
  meta <- data.frame(sample_data(ps_fam), stringsAsFactors = FALSE)
  n_g1 <- sum(meta[[group_var]] == group1, na.rm = TRUE)
  n_g2 <- sum(meta[[group_var]] == group2, na.rm = TRUE)

  cat("  ", group1, ":", n_g1, "\n", sep = "")
  cat("  ", group2, ":", n_g2, "\n", sep = "")

  if (n_g1 < 3 || n_g2 < 3) {
    cat("  WARNING: Insufficient samples per group (need >= 3)\n")
    return(NULL)
  }

  # 1. ASV overlap
  overlap_result <- calculate_asv_overlap_jaccard(ps_fam, group_var, group1, group2)

  cat("  Shared ASVs:", overlap_result$n_shared, "\n")
  cat("  ", group1, "-only:", overlap_result$n_unique_g1, "\n", sep = "")
  cat("  ", group2, "-only:", overlap_result$n_unique_g2, "\n", sep = "")
  cat("  Jaccard similarity:", round(overlap_result$jaccard, 3), "\n")

  # 2. Abundance testing for shared ASVs
  abundance_tests <- NULL
  n_sig <- 0
  pct_sig <- 0

  if (overlap_result$n_shared > 0) {
    abundance_tests <- test_shared_asv_abundances(
      ps_fam, overlap_result$shared_asvs, group_var, group1, group2
    )
    n_sig <- sum(abundance_tests$Significant, na.rm = TRUE)
    pct_sig <- 100 * n_sig / nrow(abundance_tests)
    cat("  Significant abundance differences:", n_sig, "/",
      nrow(abundance_tests), "(", round(pct_sig, 1), "%)\n",
      sep = ""
    )
  }

  # 3. Bray-Curtis partitioning
  partition_result <- partition_bray_curtis_turnover_abundance(
    ps_fam, group_var, group1, group2
  )

  cat("  BC total:", round(partition_result$bc_total, 3), "\n")
  cat("  Turnover:", round(100 * partition_result$prop_turnover, 1), "%\n", sep = "")
  cat("  Abundance:", round(100 * partition_result$prop_abundance, 1), "%\n", sep = "")

  # Determine dominant mechanism
  dominant <- if (partition_result$prop_turnover > 0.6) {
    "Turnover"
  } else if (partition_result$prop_abundance > 0.6) {
    "Abundance"
  } else {
    "Mixed"
  }

  cat("  Dominant mechanism:", dominant, "\n")

  return(list(
    family = family_name,
    n_asvs = family_subset$n_asvs,
    n_samples = nsamples(ps_fam),
    n_g1 = n_g1,
    n_g2 = n_g2,
    overlap = overlap_result,
    abundance_tests = abundance_tests,
    n_sig_abundance = n_sig,
    pct_sig_abundance = pct_sig,
    partition = partition_result,
    dominant_mechanism = dominant
  ))
}


#' Run drill-down analysis for multiple families from SIMPER results
#'
#' @param ps phyloseq object
#' @param simper_table Data frame with Family and contr_pct columns
#' @param group_var Character, grouping variable
#' @param group1 Character, first group name
#' @param group2 Character, second group name
#' @param primerset_name Character, primer set label
#' @param contr_threshold Numeric, minimum % contribution (default 5)
#' @param max_families Integer, maximum families to analyze (default 10)
#' @return Named list of drill-down results by family
run_multi_family_drilldown <- function(ps, simper_table, group_var, group1, group2,
                                       primerset_name, contr_threshold = 5,
                                       max_families = 10) {
  cat("\n=== ", primerset_name, ": MULTI-FAMILY DRILL-DOWN ===\n", sep = "")
  cat("Threshold: >", contr_threshold, "% contribution\n", sep = "")

  # Filter to qualifying families
  qualifying_families <- simper_table %>%
    filter(contr_pct >= contr_threshold) %>%
    arrange(desc(contr_pct)) %>%
    slice_head(n = max_families)

  if (nrow(qualifying_families) == 0) {
    cat("No families meet the contribution threshold\n")
    return(list())
  }

  cat("Analyzing", nrow(qualifying_families), "families\n")

  # Run drill-down for each family
  results <- list()

  for (i in 1:nrow(qualifying_families)) {
    family_name <- as.character(qualifying_families$Family[i])
    contr <- as.numeric(qualifying_families$contr_pct[i])

    cat("\n[", i, "/", nrow(qualifying_families), "] ", family_name,
      " (", round(contr, 1), "% contribution)\n",
      sep = ""
    )

    family_result <- run_family_drilldown(
      ps, family_name, group_var, group1, group2, primerset_name
    )

    if (!is.null(family_result)) {
      family_result$simper_contribution <- contr
      results[[family_name]] <- family_result
    }
  }

  cat("\n=== SUMMARY ===\n")
  cat(
    "Successfully analyzed", length(results), "/",
    nrow(qualifying_families), "families\n"
  )

  return(results)
}


#' Summarize drill-down results into publication-ready table
#'
#' @param drilldown_results List of drill-down results from run_multi_family_drilldown
#' @return Data frame with summary statistics by family
summarize_drilldown_results <- function(drilldown_results) {
  if (length(drilldown_results) == 0) {
    cat("No results to summarize\n")
    return(data.frame())
  }

  summary_list <- list()

  for (family_name in names(drilldown_results)) {
    res <- drilldown_results[[family_name]]

    summary_list[[family_name]] <- data.frame(
      Family = as.character(family_name),
      SIMPER_Contr_Pct = as.numeric(res$simper_contribution),
      N_ASVs = as.integer(res$n_asvs),
      N_Shared = as.integer(res$overlap$n_shared),
      N_Unique_G1 = as.integer(res$overlap$n_unique_g1),
      N_Unique_G2 = as.integer(res$overlap$n_unique_g2),
      Jaccard = as.numeric(res$overlap$jaccard),
      BC_Total = as.numeric(res$partition$bc_total),
      BC_Turnover_Pct = as.numeric(100 * res$partition$prop_turnover),
      BC_Abundance_Pct = as.numeric(100 * res$partition$prop_abundance),
      Dominant_Mechanism = as.character(res$dominant_mechanism),
      Pct_Sig_Abundance = as.numeric(res$pct_sig_abundance),
      stringsAsFactors = FALSE
    )
  }

  # Check if we have any results
  if (length(summary_list) == 0) {
    cat("No families had sufficient data for analysis\n")
    return(data.frame())
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  # Sort by SIMPER contribution
  summary_df <- summary_df[order(-summary_df$SIMPER_Contr_Pct), ]

  return(summary_df)
}


#' Create faceted visualization of turnover vs abundance by family
#'
#' @param summary_table Data frame from summarize_drilldown_results
#' @param title Character, plot title
#' @return ggplot object
plot_drilldown_summary <- function(summary_table, title = "") {
  if (nrow(summary_table) == 0) {
    cat("No data to plot\n")
    return(NULL)
  }

  # Reshape for stacked bars
  plot_data <- summary_table %>%
    select(Family, BC_Turnover_Pct, BC_Abundance_Pct, Dominant_Mechanism) %>%
    pivot_longer(
      cols = c(BC_Turnover_Pct, BC_Abundance_Pct),
      names_to = "Component", values_to = "Percentage"
    ) %>%
    mutate(
      Component = ifelse(Component == "BC_Turnover_Pct", "Turnover", "Abundance"),
      Family = factor(Family, levels = rev(summary_table$Family)) # Order by contribution
    )

  # Create plot
  p <- ggplot(plot_data, aes(x = Percentage, y = Family, fill = Component)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Turnover" = "#0072B2", "Abundance" = "#D55E00")) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = "Partitioning of Bray-Curtis dissimilarity",
      x = "Percentage of BC Dissimilarity",
      y = "Prey Family",
      fill = "Component"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "gray50", alpha = 0.5)

  return(p)
}

# ==============================================================================
# ASV-LEVEL TWO-PANEL VISUALIZATION FUNCTIONS
# ==============================================================================

#' Calculate ASV Category Proportions
#'
#' Calculate total reads from Shared, G1-only, and G2-only ASVs for each group
#'
#' @param ps phyloseq object
#' @param overlap_result Result from calculate_asv_overlap_jaccard()
#' @param group_var Grouping variable name in sample_data
#' @param group1 First group level
#' @param group2 Second group level
#' @return Data frame with columns Group, Category, Reads, Proportion
calculate_asv_category_proportions <- function(ps, overlap_result, group_var, group1, group2) {
  # Extract OTU matrix (ensure taxa as rows)
  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }

  # Get sample IDs for each group
  meta <- as.data.frame(sample_data(ps))
  samples_g1 <- rownames(meta)[meta[[group_var]] == group1]
  samples_g2 <- rownames(meta)[meta[[group_var]] == group2]

  # Initialize results
  results <- data.frame()

  # Process Group 1
  if (length(samples_g1) > 0) {
    # Sum reads for each category
    shared_reads <- sum(otu_mat[overlap_result$shared_asvs, samples_g1, drop = FALSE])
    g1_only_reads <- sum(otu_mat[overlap_result$unique_g1, samples_g1, drop = FALSE])
    g2_only_reads <- sum(otu_mat[overlap_result$unique_g2, samples_g1, drop = FALSE])

    total_g1 <- shared_reads + g1_only_reads + g2_only_reads

    results <- rbind(results, data.frame(
      Group = rep(group1, 3),
      Category = c("Shared", paste0(group1, "-only"), paste0(group2, "-only")),
      Reads = c(shared_reads, g1_only_reads, g2_only_reads),
      Proportion = c(shared_reads, g1_only_reads, g2_only_reads) / total_g1 * 100
    ))
  }

  # Process Group 2
  if (length(samples_g2) > 0) {
    shared_reads <- sum(otu_mat[overlap_result$shared_asvs, samples_g2, drop = FALSE])
    g1_only_reads <- sum(otu_mat[overlap_result$unique_g1, samples_g2, drop = FALSE])
    g2_only_reads <- sum(otu_mat[overlap_result$unique_g2, samples_g2, drop = FALSE])

    total_g2 <- shared_reads + g1_only_reads + g2_only_reads

    results <- rbind(results, data.frame(
      Group = rep(group2, 3),
      Category = c("Shared", paste0(group1, "-only"), paste0(group2, "-only")),
      Reads = c(shared_reads, g1_only_reads, g2_only_reads),
      Proportion = c(shared_reads, g1_only_reads, g2_only_reads) / total_g2 * 100
    ))
  }

  # Set factor levels for Category (for stacking order: G2-only, Shared, G1-only)
  results$Category <- factor(results$Category,
    levels = c(paste0(group2, "-only"), "Shared", paste0(group1, "-only"))
  )
  results$Group <- factor(results$Group, levels = c(group1, group2))

  return(results)
}

#' Extract Top SIMPER ASVs with Taxonomy and Category
#'
#' Extract top N ASVs from SIMPER with taxonomy labels and category assignment
#'
#' @param simper_result SIMPER result object from vegan::simper()
#' @param comparison_name Name of comparison in SIMPER result (e.g., "Male_Female")
#' @param ps phyloseq object for taxonomy lookup
#' @param overlap_result Result from calculate_asv_overlap_jaccard()
#' @param group1 First group name
#' @param group2 Second group name
#' @param top_n Number of top ASVs to extract (default 20)
#' @return Data frame with ASV, Genus, Family, Contribution_Pct, Category, Display_Label
extract_top_simper_asvs <- function(simper_result, comparison_name, ps, overlap_result,
                                    group1, group2, top_n = 20) {
  # Extract SIMPER data for the comparison
  simper_comp <- simper_result[[comparison_name]]

  # Get ASV names, contributions, and group averages
  asv_names <- simper_comp$species
  avg_contrib <- simper_comp$average
  overall_contrib <- simper_comp$overall
  ava <- simper_comp$ava # average abundance in group A
  avb <- simper_comp$avb # average abundance in group B

  # Calculate contribution percentage
  contrib_pct <- (avg_contrib / overall_contrib) * 100

  # Create initial data frame
  simper_df <- data.frame(
    ASV = asv_names,
    Contribution = avg_contrib,
    Contribution_Pct = contrib_pct,
    Avg_G1 = ava,
    Avg_G2 = avb,
    stringsAsFactors = FALSE
  )

  # Sort by contribution and take top N
  simper_df <- simper_df[order(-simper_df$Contribution_Pct), ]
  simper_df <- head(simper_df, top_n)

  # Get taxonomy for these ASVs
  tax <- as.data.frame(tax_table(ps))
  simper_df$Genus <- tax[simper_df$ASV, "Genus"]
  simper_df$Family <- tax[simper_df$ASV, "Family"]

  # Clean QIIME2 prefixes
  simper_df$Genus <- sub("^g__", "", simper_df$Genus)
  simper_df$Family <- sub("^f__", "", simper_df$Family)

  # Replace empty/NA genus with truncated ASV ID
  simper_df$Genus <- ifelse(
    is.na(simper_df$Genus) | simper_df$Genus == "" | simper_df$Genus == "uncultured",
    paste0(substr(simper_df$ASV, 1, 12), "..."),
    simper_df$Genus
  )

  # Categorize each ASV
  simper_df$Category <- sapply(simper_df$ASV, function(asv) {
    if (asv %in% overlap_result$shared_asvs) {
      return("Shared")
    } else if (asv %in% overlap_result$unique_g1) {
      return(paste0(group1, "-only"))
    } else if (asv %in% overlap_result$unique_g2) {
      return(paste0(group2, "-only"))
    } else {
      return("Unknown")
    }
  })

  # Create display labels: "Genus (Family)"
  simper_df$Display_Label <- paste0(simper_df$Genus, " (", simper_df$Family, ")")

  # Make labels unique by appending ASV suffix when duplicates exist
  label_counts <- table(simper_df$Display_Label)
  duplicated_labels <- names(label_counts)[label_counts > 1]

  if (length(duplicated_labels) > 0) {
    for (dup_label in duplicated_labels) {
      dup_indices <- which(simper_df$Display_Label == dup_label)
      # Append short ASV identifier for duplicates
      for (i in seq_along(dup_indices)) {
        idx <- dup_indices[i]
        asv_suffix <- substr(simper_df$ASV[idx], nchar(simper_df$ASV[idx]) - 5, nchar(simper_df$ASV[idx]))
        simper_df$Display_Label[idx] <- paste0(simper_df$Display_Label[idx], " [", asv_suffix, "]")
      }
    }
  }

  # Set factor levels for Category (same as proportion function)
  simper_df$Category <- factor(simper_df$Category,
    levels = c(paste0(group2, "-only"), "Shared", paste0(group1, "-only"))
  )

  return(simper_df)
}

#' Plot ASV Category Stacked Bars
#'
#' Create Panel A: stacked bar chart of read proportions by ASV category
#'
#' @param proportion_data Data frame from calculate_asv_category_proportions()
#' @param group1 First group name
#' @param group2 Second group name
#' @param title Plot title (default "")
#' @return ggplot object
plot_asv_category_bars <- function(proportion_data, group1, group2, title = "") {
  # Define color palette
  color_palette <- c(
    "Shared" = "gray60",
    setNames("#0072B2", paste0(group1, "-only")), # Blue for group1
    setNames("#D55E00", paste0(group2, "-only")) # Red/orange for group2
  )

  p <- ggplot(proportion_data, aes(x = Group, y = Proportion, fill = Category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Proportion)),
      position = position_stack(vjust = 0.5),
      size = 3.5, color = "white", fontface = "bold"
    ) +
    scale_fill_manual(values = color_palette, name = "ASV Category") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    labs(
      title = ifelse(title != "", title, "ASV Category Composition"),
      subtitle = "Proportion of reads from shared and unique ASVs",
      x = NULL,
      y = "Proportion of Reads (%)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  return(p)
}

#' Plot Top ASV Drivers
#'
#' Create Panel B: horizontal bar plot of top ASVs driving dissimilarity
#'
#' @param top_asv_data Data frame from extract_top_simper_asvs()
#' @param group1 First group name
#' @param group2 Second group name
#' @param title Plot title (default "")
#' @return ggplot object
plot_top_asv_drivers <- function(top_asv_data, group1, group2, title = "") {
  # Define same color palette as Panel A
  color_palette <- c(
    "Shared" = "gray60",
    setNames("#0072B2", paste0(group1, "-only")),
    setNames("#D55E00", paste0(group2, "-only"))
  )

  # Reverse order for plotting (top contributor at top)
  top_asv_data$Display_Label <- factor(top_asv_data$Display_Label,
    levels = rev(top_asv_data$Display_Label)
  )

  p <- ggplot(top_asv_data, aes(x = Contribution_Pct, y = Display_Label, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_palette, name = "ASV Category") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      title = ifelse(title != "", title, "Top ASVs Driving Dissimilarity"),
      subtitle = "From SIMPER analysis (Bray-Curtis)",
      x = "Contribution to Dissimilarity (%)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "none", # Legend inherited from Panel A
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  return(p)
}

#' Create ASV Two-Panel Figure
#'
#' Wrapper function orchestrating complete two-panel ASV visualization
#'
#' @param ps phyloseq object
#' @param simper_result SIMPER result object
#' @param comparison_name Name of comparison in SIMPER result
#' @param group_var Grouping variable name
#' @param group1 First group level
#' @param group2 Second group level
#' @param main_title Overall figure title (default "")
#' @param top_n Number of top ASVs to show (default 15)
#' @return patchwork combined figure
create_asv_twopanel_figure <- function(ps, simper_result, comparison_name,
                                       group_var, group1, group2,
                                       main_title = "", top_n = 15) {
  # Calculate ASV overlap
  cat("Calculating ASV overlap between", group1, "and", group2, "...\n")
  overlap_result <- calculate_asv_overlap_jaccard(ps, group_var, group1, group2)

  # Print summary
  cat("\nASV Overlap Summary:\n")
  cat("  Shared:", length(overlap_result$shared_asvs), "\n")
  cat(" ", group1, "-only:", length(overlap_result$unique_g1), "\n")
  cat(" ", group2, "-only:", length(overlap_result$unique_g2), "\n")
  cat("  Jaccard similarity:", round(overlap_result$jaccard, 3), "\n\n")

  # Calculate category proportions
  proportion_data <- calculate_asv_category_proportions(
    ps, overlap_result,
    group_var, group1, group2
  )

  # Extract top SIMPER ASVs
  top_asv_data <- extract_top_simper_asvs(
    simper_result, comparison_name, ps,
    overlap_result, group1, group2, top_n
  )

  cat("Extracted top", nrow(top_asv_data), "ASVs from SIMPER\n")

  # Create panels
  panel_a <- plot_asv_category_bars(proportion_data, group1, group2,
    title = "A. ASV Category Composition"
  )

  panel_b <- plot_top_asv_drivers(top_asv_data, group1, group2,
    title = paste0(
      "B. Top ", nrow(top_asv_data),
      " ASVs Driving Dissimilarity"
    )
  )

  # Combine with patchwork
  combined <- panel_a + panel_b +
    plot_layout(ncol = 2, widths = c(1, 2), guides = "collect") +
    plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    ) &
    theme(legend.position = "bottom")

  return(combined)
}
