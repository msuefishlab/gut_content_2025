# Mormyrid Gut Content Analysis: Revision Plan

## Project Overview

This document guides revision of R markdown analysis files for a DNA metabarcoding study of mormyrid electric fish gut contents in Gabon. The study uses two independent COI primer sets (PS1: mlCOIintF/BR2; PS2: BF2/BR1) for cross-validation.

**Core biological question:** Does EOD waveform polymorphism (biphasic vs. triphasic) in *Paramormyrops kingsleyae* associate with dietary niche partitioning?

**Key design consideration:** BP and TP *P. kingsleyae* are never sympatric - waveform is confounded with geography. The analysis must separate these effects.

---

## Analysis Framework: Four Biological Levels

Structure the R markdown to flow through these levels, each with its own section.

---

### LEVEL 1: Mormyrid Community (Balé Creek)
*Scope: 9 species, 35 individuals from single site - controls for geography*

#### Existing Analyses (keep, may need restructuring):
- [x] Alpha diversity among species (Shannon, ANOVA)
- [x] Beta diversity among species (Bray-Curtis, PERMANOVA)
- [x] Diet overlap indices (Schoener, Pianka) - calculated but not statistically tested
- [x] Pairwise species comparisons (pairwiseAdonis)
- [x] PCoA visualization

#### NEW Analyses Needed:

**1.1 Prey Composition Summary**
```r
# Aggregate ASVs to Family or Order level
# Calculate relative abundance per species
# Create stacked bar chart or heatmap showing:
#   - Rows: Species (or individual fish)
#   - Columns/Colors: Prey families
# Report top 5-10 dominant prey families with percentages
```

**1.2 Phylogenetic Correlation Test**
```r
# Need: Phylogenetic distance matrix for the 9 species
# (derive from cytb tree in Arnegard et al. 2010 or similar)
# 
# Mantel test:
# mantel(diet_bray_curtis, phylo_distance, method = "pearson", permutations = 999)
#
# Tests: Do more closely related species have more similar diets?
```

**1.3 Within-Paramormyrops vs. Between-Genus Overlap Test**
```r
# From existing overlap matrices (Schoener or Pianka):
# 
# 1. Extract pairwise overlaps for:
#    - Paramormyrops-Paramormyrops pairs
#    - Paramormyrops-other genus pairs (or other-other pairs)
#
# 2. Compare distributions:
#    wilcox.test(within_paramormyrops_overlap, between_genus_overlap)
#    # or t.test() if normally distributed
#
# Tests: Do Paramormyrops species overlap more with each other?
# (relevant to competitive exclusion)
```

**1.4 SIMPER Analysis - Species**
```r
# library(vegan)
# simper_species <- simper(asv_table, species_groups, permutations = 999)
# summary(simper_species)
#
# Identifies which prey taxa contribute most to between-species differences
```

---

### LEVEL 2: BP vs. TP Globally (All Mormyrids)
*Scope: All BP (n≈63-65) vs. all TP (n≈34) fish across all sites/species*

#### Existing Analyses (keep):
- [x] Alpha diversity BP vs. TP (Kruskal-Wallis or t-test)
- [x] Beta diversity BP vs. TP (PERMANOVA) - R² ≈ 0.035

#### NEW Analyses Needed:

**2.1 Prey Composition by Waveform**
```r
# Aggregate to Family/Order level
# Create visualization comparing BP vs. TP:
#   - Side-by-side stacked bars, OR
#   - Heatmap with BP and TP columns
# 
# Calculate mean ± SE relative abundance of top taxa for each waveform
```

**2.2 SIMPER Analysis - Waveform**
```r
# simper_waveform <- simper(asv_table, waveform_type, permutations = 999)
# summary(simper_waveform)
#
# Which prey families drive BP vs. TP compositional differences?
```

**2.3 Explicit Alpha Diversity Comparison**
```r
# The existing analysis shows TP has higher diversity - make this explicit
# Report: mean Shannon ± SD for BP vs. TP
# Effect size (Cohen's d or similar)
```

---

### LEVEL 3: *Paramormyrops kingsleyae* - General Factors
*Scope: ~47 P. kingsleyae from 4 rivers (Bambomo, Mengono = BP; Apassa, Iboundji = TP)*

#### Existing Analyses (keep):
- [x] Sex effects on Shannon diversity (ANCOVA with body size)
- [x] Body size effects on Shannon diversity (linear regression)

#### NEW Analyses Needed:

**3.1 Prey Composition by Population**
```r
# Stacked bar chart with 4 populations (rivers) side by side
# Group visually by region (North: Mengono + Iboundji; South: Apassa + Bambomo)
# OR group by waveform (BP: Bambomo + Mengono; TP: Apassa + Iboundji)
```

**3.2 Parasitism Effects**
```r
# IF parasitism data exists in metadata:
#
# Alpha diversity:
# kruskal.test(shannon ~ parasitism_status, data = pk_data)
# # or ANCOVA if continuous parasite load
#
# Beta diversity:
# adonis2(asv_table ~ parasitism_status, data = pk_data)
#
# IF NO parasitism data: Note this as a limitation/future direction
```

---

### LEVEL 4: BP vs. TP Within *P. kingsleyae*
*The core question - tests dietary niche partitioning hypothesis*

#### Existing Analyses (keep, improve presentation):
- [x] Two-way ANOVA: Shannon ~ Region + Waveform (Region sig, Waveform NS)
- [x] Two-way PERMANOVA: Bray-Curtis ~ Region + Waveform (Both sig, waveform R² ≈ 0.10)
- [x] Prey selectivity (Ivlev's E) - strong non-random feeding
- [x] Selectivity profiles BP vs. TP (PERMANOVA on E profiles - NS)
- [x] PCoA visualization by region and waveform

#### NEW Analyses Needed:

**4.1 Prey Composition BP vs. TP in P. kingsleyae**
```r
# Create figure with 4 panels or facets:
#   - Bambomo (BP, South)
#   - Apassa (TP, South)  
#   - Mengono (BP, North)
#   - Iboundji (TP, North)
#
# Arrange to facilitate BP vs. TP comparison within each region
```

**4.2 Waveform Similarity Across Sites**
```r
# Extract pairwise Bray-Curtis distances
# Categorize each pair:
#   - "Same waveform, same region" (e.g., within Bambomo)
#   - "Same waveform, different region" (e.g., Bambomo vs. Mengono)
#   - "Different waveform, same region" (e.g., Bambomo vs. Apassa)
#   - "Different waveform, different region" (e.g., Bambomo vs. Iboundji)
#
# Compare mean distances across categories (ANOVA or Kruskal-Wallis)
# 
# Key test: Is same-waveform-different-region distance < different-waveform-same-region?
# If yes: waveform predicts diet better than geography
```

**4.3 Variance Partitioning**
```r
# Use varpart() from vegan, or compare sequential vs. marginal R² from PERMANOVA
#
# adonis2(asv ~ Region + Waveform, data = pk_data, by = "margin")
# 
# Report:
#   - Variance explained by Region alone
#   - Variance explained by Waveform alone  
#   - Shared variance
#   - Residual
```

**4.4 Mean Selectivity Comparison**
```r
# For each fish, calculate mean |E| (absolute Ivlev's electivity)
# This measures how "choosy" each fish is overall
#
# pk_data$mean_abs_E <- rowMeans(abs(ivlev_matrix), na.rm = TRUE)
#
# Compare BP vs. TP:
# wilcox.test(mean_abs_E ~ waveform, data = pk_data)
#
# Are BP or TP more selective overall?
```

**4.5 Integration Point (interpretive, not new stats)**
```r
# In the Results narrative, explicitly connect:
#   - Beta diversity: BP and TP differ in WHAT they eat
#   - Selectivity profiles: BP and TP do NOT differ in HOW they select
#   - Conclusion: Dietary differences likely reflect prey availability, not foraging strategy
```

---

## Cross-Validation Approach

For ALL analyses above, run on both primer sets independently. Report:
- Whether conclusions agree (direction and significance)
- Correlation coefficients for continuous metrics
- Mantel test r-values for distance matrices

Move detailed cross-validation metrics to supplementary material or a summary table. Don't lead sections with validation - lead with biology.

---

## Output Structure Suggestion

```
mormyrid_gut_analysis_revised.Rmd
│
├── 0. Setup & Data Loading
│   ├── Load packages
│   ├── Import QIIME2 artifacts (both primer sets)
│   ├── Filter host sequences
│   ├── Rarefy to 5000 reads
│   └── Create phyloseq objects
│
├── 1. Mormyrid Community at Balé Creek
│   ├── 1.1 Prey composition (NEW)
│   ├── 1.2 Alpha diversity among species
│   ├── 1.3 Beta diversity among species
│   ├── 1.4 Diet overlap indices + statistical test (MODIFIED)
│   ├── 1.5 Phylogenetic correlation (NEW)
│   └── 1.6 SIMPER - species drivers (NEW)
│
├── 2. BP vs. TP Mormyrids Globally
│   ├── 2.1 Prey composition by waveform (NEW)
│   ├── 2.2 Alpha diversity comparison
│   ├── 2.3 Beta diversity comparison
│   └── 2.4 SIMPER - waveform drivers (NEW)
│
├── 3. P. kingsleyae - General Factors
│   ├── 3.1 Prey composition by population (NEW)
│   ├── 3.2 Sex effects
│   ├── 3.3 Body size effects
│   └── 3.4 Parasitism effects (NEW, if data available)
│
├── 4. BP vs. TP Within P. kingsleyae
│   ├── 4.1 Prey composition (NEW)
│   ├── 4.2 Alpha diversity (Region + Waveform)
│   ├── 4.3 Beta diversity (Region + Waveform)
│   ├── 4.4 Variance partitioning (NEW)
│   ├── 4.5 Waveform similarity across sites (NEW)
│   ├── 4.6 Prey selectivity analysis
│   ├── 4.7 Mean selectivity by waveform (NEW)
│   └── 4.8 Synthesis: composition vs. selectivity
│
├── 5. Cross-Validation Summary
│   └── Table of concordance across all tests
│
└── 6. Session Info
```

---

## Data Requirements Checklist

Before running revisions, confirm availability of:

- [ ] Phylogenetic tree or distance matrix for 9 Balé Creek species
- [ ] Parasitism status in specimen metadata
- [ ] Environmental sampling data (for selectivity analysis - appears to exist)
- [ ] Both primer set phyloseq objects with matching sample IDs

---

## Notes for Claude Code

1. **Preserve existing working code** - refactor into this structure, don't rewrite from scratch
2. **Run both primer sets** - every analysis should have PS1 and PS2 versions
3. **Consistent visualization style** - use ggplot2 with consistent color palettes for waveform (BP/TP) and region (North/South)
4. **Table outputs** - create summary tables suitable for manuscript (kable or gt format)
5. **Commenting** - add comments explaining the biological question each analysis addresses
