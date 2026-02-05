# Filtering Loss Analysis

## Overview

This analysis quantifies the loss of ASVs and reads during Family-level taxonomic filtering in the gut content metabarcoding pipeline. The goal is to determine whether differential filtering rates between primer sets could explain observed compositional differences.

## Background

The pipeline filters ASVs to retain only those assigned to at least **Family level** (both Order AND Family must be non-NA). This strict filtering criterion (implemented in `code/06_Generate_Output/filterVsearch_nbClassifier.R`, lines 48-50) means that "100% family assignment" in downstream analyses is by construction—unassigned ASVs were removed upstream.

This analysis quantifies:
1. **What was lost** during filtering (ASV counts and read abundance)
2. **Whether PS1 and PS2 were affected differently** (differential filtering)
3. **Why ASVs were removed** (taxonomic assignment level)
4. **Which Orders are most affected** by the Family-level requirement

## Input Files

### Before Filtering (raw clustered data)
- `output_data/03_Clustered_Data/primerset1_all_p985_table.qza`
- `output_data/03_Clustered_Data/primerset2_all_p985_table.qza`

### After Filtering (family-assigned only)
- `output_data/06_Generate_Output/primerset1_all_p985_table_filtd.qza`
- `output_data/06_Generate_Output/primerset2_all_p985_table_filtd.qza`

### Taxonomy Assignments
- `output_data/06_Generate_Output/primerset1_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv`
- `output_data/06_Generate_Output/primerset2_all_p985_taxa_VsearchOnly_p95_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv`

## Execution

### On HPCC (with Singularity)

```bash
# Ensure gut_contents.env is sourced
source gut_contents.env

# Run analysis
bash code/07_Analysis/run_filtering_loss_analysis.sh
```

### Locally (without Singularity)

```bash
# Run analysis (sets root automatically)
bash code/07_Analysis/run_filtering_loss_analysis_local.sh
```

**Note:** Local execution requires R packages (tidyverse, qiime2R, patchwork) installed in your local R environment.

### Expected Runtime
2-3 minutes (loading .qza files, parsing taxonomy, generating plots)

## Output Files

All outputs are saved to: `output_data/07_Analysis/filtering_loss_analysis/`

### Tables (CSV format)
1. **`global_filtering_summary.csv`** - Primer-level totals (ASVs and reads before/after/lost)
2. **`per_sample_filtering_loss.csv`** - Sample-level detail (each sample's loss %)
3. **`filtered_asv_taxonomy_breakdown.csv`** - Removed ASVs by category
4. **`lost_reads_by_order.csv`** - Orders most affected by Family-level requirement
5. **`statistical_tests_summary.csv`** - All test results with p-values
6. **`removed_asvs_full_list.csv`** - Complete list of removed ASVs with taxonomy

### Figures (PDF format)
1. **`filtering_loss_comparison.pdf`** - PS1 vs PS2 global loss (bar plot)
2. **`per_sample_loss_distribution.pdf`** - Boxplot + histogram of sample-level loss
3. **`lost_taxonomy_breakdown.pdf`** - Stacked bar of removal categories
4. **`lost_reads_by_order.pdf`** - Top Orders with Order-only assignments
5. **`removed_asv_abundance_distribution.pdf`** - Abundance distribution of lost ASVs

### Report
**`report/filtering_loss_report.md`** - Comprehensive markdown report with findings

## Analysis Components

### 1. Global Filtering Loss
Calculates total ASVs and reads lost for each primer set, both as counts and percentages.

### 2. Per-Sample Filtering Loss
Examines variation in filtering loss across individual samples to identify outliers.

### 3. Taxonomy Breakdown
Categorizes removed ASVs by taxonomic assignment level:
- **Unassigned** - No assignment at any level
- **Order-only (no Family)** - Assigned to Order but not Family
- **Class-only (no Order)** - Assigned to Class but not Order
- **Phylum-only** - Assigned to Phylum but not Class
- **Kingdom-only** - Assigned to Kingdom but not Phylum
- **Other** - Edge cases

### 4. Order-Level Loss
Identifies which Orders are most affected by the Family-level filtering requirement (i.e., Orders with many ASVs assigned at Order-level but lacking Family assignments).

### 5. Statistical Tests
- **Two-proportion z-test**: Compares global % reads lost between PS1 and PS2
- **Wilcoxon rank-sum test**: Compares per-sample loss distributions
- **Chi-square test**: Tests whether category distributions differ by primer set

## Interpretation Guidelines

### If PS1 and PS2 have similar filtering loss rates (difference < 5%)
**Interpretation:** Compositional differences are likely **real biological differences** in primer amplification specificity, not taxonomic assignment artifacts.

### If PS1 has substantially higher loss than PS2 (difference > 10%)
**Interpretation:** The observed compositional patterns may be **partially driven by differential filtering**—some taxa may have been disproportionately removed due to poorer Family-level assignments in one primer set.

### If Order-only assignments are concentrated in specific Orders
**Interpretation:** These Orders represent taxonomic groups where the reference database lacks Family-level resolution. These are candidates for:
- Manual BLAST verification
- Reference database improvement
- Relaxing filtering criteria for specific taxa

## Integration with Pipeline

This analysis is part of the broader investigation into primer set differences:

1. **assignment_diagnostics_analysis.R** - Analyzes RETAINED ASVs
2. **filtering_loss_analysis.R** (this script) - Analyzes REMOVED ASVs
3. **local_barcode_contribution_analysis.R** - Investigates locally-barcoded ASV patterns
4. Main diversity analyses - Interpret diversity metrics in context of filtering loss

## Troubleshooting

### Error: "Input file not found"
Check that all input files exist. The analysis requires:
- Both primer sets' before/after feature tables (.qza files)
- Both primer sets' taxonomy files (.tsv files)

### Error: "'root' environment variable not set"
Source the environment file first:
```bash
source gut_contents.env
```

### Warning: "Some samples lost >50% of reads"
This indicates samples with extreme filtering loss. Review per-sample results to identify which samples are affected and why.

### Warning: "PS1 and PS2 show >10% difference in filtering loss rate"
This suggests differential filtering that could partially explain compositional differences. Review the statistical tests and category breakdown to understand the source of the difference.

## Dependencies

- R packages: tidyverse, qiime2R, patchwork, knitr
- Singularity container: `${rimage}` (defined in gut_contents.env)
- Environment file: `gut_contents.env`

## Author

Created: 2026-02-04
