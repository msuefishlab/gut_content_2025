#===== DATA LOADING =====
  
  #' Load primer set data into phyloseq object
  #'
  #' @param primerset_num Integer (1 or 2) indicating which primer set
  #' @param use_no_host Logical, use NO_HOST filtered data (default TRUE)
  #' @return phyloseq object
  load_primerset_data <- function(primerset_num, use_no_host = TRUE) {
    ps_name <- paste0("primerset", primerset_num)
    
    # Construct file paths
    metadata_file <- file.path(root, "input_data", "metadata",
                               "wholestudy_metadata_cleaned_jan_2026.tsv")
    
    if (use_no_host) {
      table_file <- file.path(root, "output_data", "06_Generate_Output",
                              paste0(ps_name, "_all_p985_table_filtd_NO_HOST.rarefied.qza"))
    } else {
      table_file <- file.path(root, "output_data", "06_Generate_Output",
                              paste0(ps_name, "_all_p985_table_filtd.rarefied.qza"))
    }
    
    taxa_file <- file.path(root, "output_data", "06_Generate_Output",
                           paste0(ps_name, "_all_p985_taxa_filtd_ALL.qza"))
    
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
  div[["Sample_Name"]]<-meta$Sample_Name
  #div$SampleID <- rownames(div)
  
  # Test normality
  shapiro_p <- shapiro.test(div$Shannon)$p.value
  
  # Choose appropriate test
  if (shapiro_p < 0.05) {
    # Non-parametric
    test_result <- kruskal.test(as.formula(paste("Shannon ~", group_var)),
                                data = div)
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
                         permutations = 999)
  
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
  p <- ggplot(div_combined, aes(x = .data[[group_var]], y = Shannon,
                                fill = PrimerSet)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.7) +
    scale_fill_manual(values = c("PS1" = ps1_color, "PS2" = ps2_color)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = group_var, y = "Shannon Diversity Index",
         fill = "Primer Set")
  
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
    p <- ggplot(pcoa_combined, aes(x = PCoA1, y = PCoA2,
                                   color = .data[[color_var]])) +
      geom_point(size = 3, alpha = 0.7)
  } else {
    p <- ggplot(pcoa_combined, aes(x = PCoA1, y = PCoA2,
                                   color = .data[[color_var]],
                                   shape = .data[[shape_var]])) +
      geom_point(size = 3, alpha = 0.7)
  }
  
  p <- p +
    facet_wrap(~ PrimerSet, scales = "free") +
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
    return(list(cor_test = NULL, n_shared = length(shared),
                message = "Insufficient shared samples"))
  }
  
  # Subset and match order
  div1_shared <- div1[match(shared, div1$Sample_Name), ]
  div2_shared <- div2[match(shared, div2$Sample_Name), ]
  
  # Pearson correlation
  cor_test <- cor.test(div1_shared$Shannon, div2_shared$Shannon,
                       method = "pearson")
  
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
  mantel_result <- mantel(dist1, dist2, method = "pearson",
                          permutations = 999)
  
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
  
  # Get samples to keep using the rownames which match sample_names(ps_obj)
  samples_to_keep <- meta[meta$Species %in% species_to_keep,]$Sample_Name
  
  # Additional safety check
  if (length(samples_to_keep) == 0) {
    cat("\n⚠ WARNING: No samples remain after filtering!\n")
    return(NULL)
  }
  
  # Verify these sample names actually exist in the phyloseq object
  actual_samples <- ps_obj@sam_data$Sample_Name
  samples_to_keep <- intersect(samples_to_keep, actual_samples)
  
  s2k<-rownames(meta)[ps_obj@sam_data$Sample_Name %in% samples_to_keep]
  
  if (length(s2k) == 0) {
    cat("\n⚠ WARNING: Sample name mismatch - no samples to keep!\n")
    cat("This suggests metadata rownames don't match phyloseq sample names.\n")
    return(NULL)
  }
  
  # Prune the phyloseq object
  ps_filt <- prune_samples(s2k, ps_obj)
  
  # CRITICAL FIX: Explicitly rebuild sample_data with correct rownames
  # This ensures metadata rownames match sample_names after all our manipulations
  meta_filt <- data.frame(sample_data(ps_filt), stringsAsFactors = FALSE)
  rownames(meta_filt) <- sample_names(ps_filt)
  sample_data(ps_filt) <- sample_data(meta_filt)
  
  cat("Retained", length(species_to_keep), "species,",
      nsamples(ps_filt), "samples\n")
  
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
  if (r + p == 0) return(NA_real_)
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
  
  if (sum_ratio == 0) return(rep(NA_real_, length(r)))
  
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
  if (denominator == 0) return(NA_real_)
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
    scale_fill_viridis_c(option = "plasma",
                         name = "Contribution\nto Dissimilarity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = paste0("Top Prey Contributing to Species Differences (",
                       primerset_name, ")"),
         x = "Species Comparison", y = "Prey Family")

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
    scale_fill_gradient2(low = avoidance_color, mid = neutral_color,
                         high = preference_color, midpoint = 0,
                         limits = c(-1, 1), na.value = "gray90",
                         name = "Ivlev's E") +
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
    facet_wrap(~ Site, ncol = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    labs(title = "Prey Availability vs. Use by Site",
         x = "", y = "Proportion", fill = "Prey Family") +
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
      tip  # Return original if can't parse
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
    query_clean <- gsub("'", "", query)  # Remove quotes for matching

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
