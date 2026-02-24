# Mormyrid Gut Content Analysis Pipeline

Metabarcoding pipeline for analyzing COI (Cytochrome Oxidase I) sequences from fish gut contents using QIIME2 and a custom COInr reference database merged with Schmidt lab barcodes.

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Container Management](#container-management)
- [Pipeline Steps](#pipeline-steps)
- [Technical Specifications](#technical-specifications)
- [Data Organization](#data-organization)
- [Troubleshooting](#troubleshooting)

## Pipeline Overview

The analysis follows a sequential structure with numbered directories under `code/`:

1. **Filter and Denoise** - Load data, trim primers, run DADA2
2. **Merge Data** - Filter low-read samples, merge runs/primers
3. **Cluster Reads** - Cluster ASVs at 98.5% identity, collapse to fish-level
4. **Reference Database** - Build custom COInr + Schmidt database, train classifier
5. **Taxonomic Assignment** - Assign taxonomy with vsearch and naive Bayes
6. **Generate Output** - Filter, combine, host-filter, rarefy, export tables
7. **Analysis** - Diagnostics, ecological analyses, and manuscript figures

## Container Management

**IMPORTANT:** Singularity build commands must be run on HPCC, not locally.

### QIIME2 Container

```bash
cd images/
singularity build --remote qiime2.sif qiime2_mod.def
```

### R Container

```bash
docker build --platform linux/amd64 -t r4guts .
docker tag r4guts jasongallant/r4guts
docker push jasongallant/r4guts
singularity build r4guts.sif docker://jasongallant/r4guts:latest
```

### mkCOInr Container

```bash
docker build --platform linux/amd64 -t mkcoinr -f mkcoinr.Dockerfile .
docker tag mkcoinr jasongallant/mkcoinr
docker push jasongallant/mkcoinr
singularity build mkcoinr.sif docker://jasongallant/mkcoinr:latest
```

## Pipeline Steps

### 01_Filter_And_Denoise

Processes raw Illumina paired-end sequencing data through quality filtering and denoising.

**What it does:**
- Loads raw FASTQ files into QIIME2 artifacts
- Trims primer sequences from paired-end reads
- Runs DADA2 denoising algorithm to identify ASVs
- Processes both `main_run` and `repeat_run` with `primerset1` and `primerset2`

**How to run:**
```bash
bash code/01_Filter_And_Denoise/01_load_data.sh
bash code/01_Filter_And_Denoise/02_trim_primers.sh
bash code/01_Filter_And_Denoise/03_inspect_trim.sh
bash code/01_Filter_And_Denoise/04_denoise.sh
bash code/01_Filter_And_Denoise/05_inspect_denoise.sh
```

All scripts submit SLURM jobs via `.sb` files. Check logs in `output_data/slurm_logs/`.

### 02_MergeData

Filters and merges denoised data across sequencing runs and primer sets.

**What it does:**
- Filters out samples below read count thresholds (`ps1_lo_thresh=1`, `ps2_lo_thresh=10` defined in `gut_contents.env`)
- Merges feature tables and representative sequences across `main_run` and `repeat_run`
- Produces a single merged table per primer set

**How to run:**
```bash
bash code/02_MergeData/01_remove_low_reads.sh
bash code/02_MergeData/02_merge_data.sh
```

### 03_Cluster_Reads

Clusters similar ASVs to reduce dataset complexity and collapses multiple PCR replicates to the fish level.

**What it does:**
- Clusters ASVs at 98.5% identity using QIIME2 vsearch de-novo clustering
- Collapses per-sample tables into per-fish tables by summing reads across replicates (grouped by `fish_id` from sample metadata)

**How to run:**
```bash
bash code/03_Cluster_Reads/01_run_cluster.sh
bash code/03_Cluster_Reads/02_collapse_to_fish.sh
```

### 04_Reference_Database

Creates a custom reference database for taxonomic assignment and trains a naive Bayes classifier.

**Key workflow:**
1. Download COInr database from Zenodo (2022-05-06 version)
2. Filter COInr to Metazoa at family level minimum using `select_taxa.pl`
3. Process Schmidt lab's 429 vouchered invertebrate barcodes from Gabon
4. Convert Schmidt barcodes to mkCOInr TSV format and add taxonomy IDs
5. Merge and dereplicate COInr + Schmidt sequences using `pool_and_dereplicate.pl`
6. Restrict to Leray primer amplicon region (280–345 bp) using `select_region.pl`
7. Format as QIIME2 artifact using `format_db.pl`
8. Train QIIME2 naive Bayes classifier

**How to run:**
```bash
# See code/04_Reference_Database/README.md for full step-by-step database build
bash code/04_Reference_Database/01_run_classifier_training.sh
```

See `code/04_Reference_Database/README.md` for the complete database construction workflow.

Uses mkCOInr Perl scripts extensively via Singularity container.

### 05_Taxanomic_Assignment

Assigns taxonomy to clustered ASVs using two complementary approaches.

**What it does:**
- **vsearch:** Exact/near-exact sequence matching for high-confidence assignments (80% identity, 94% coverage)
- **Naive Bayes classifier:** Probabilistic assignment for sequences that vsearch cannot confidently assign

**How to run:**
```bash
bash code/05_Taxanomic_Assignment/01_run_vsearch.sh
bash code/05_Taxanomic_Assignment/02_run_bayes_classifier.sh
```

### 06_Generate_Output

Selects the best taxonomy assignment per ASV, filters features, groups by fish, removes host reads, rarefies, and exports final tables.

**What it does:**
1. Exports vsearch and naive Bayes taxonomy results to TSV
2. Runs `filterVsearch_nbClassifier.R` to select the best assignment for each ASV (vsearch preferred where confident, naive Bayes otherwise)
3. Filters QIIME2 taxonomy and sequence artifacts to the selected feature IDs
4. Merges the vsearch and naive Bayes filtered taxonomy into a single combined artifact
5. Filters the clustered feature table to retained features
6. Groups samples by `fish_id` (summing reads across PCR replicates)
7. Removes host sequences (Chordata/Vertebrata/Actinopterygii/Craniata)
8. Generates rarefaction curves (500–8000 reads) to guide depth selection
9. Rarefies the host-filtered table (default 3000 reads — adjust based on `_NO_HOST.qzv`)
10. Generates QIIME2 taxa barplots (unrarefied and rarefied)

**How to run:**
```bash
bash code/06_Generate_Output/01_run_filter_and_combine.sh
bash code/06_Generate_Output/02_generate_output_host_filtered.sh
```

**Key outputs in `output_data/06_Generate_Output/`:**
- `{primer}_all_p985_taxa_filtd_ALL.qza` — combined taxonomy artifact
- `{primer}_all_p985_table_filtd_NO_HOST.qza` — host-filtered, fish-grouped feature table
- `{primer}_all_p985_table_filtd_NO_HOST.rarefied.qza` — rarefied prey-only table
- `{primer}_taxa_bar_plots_NO_HOST_no_rarefy.qzv` / `.rarefied.qzv` — QIIME2 visualizations

### 07_Analysis

Diagnostics and ecological analyses. Scripts in this directory run locally (no HPCC required for most), using system R or the R container.

#### 01_Diagnostics

Quality-check analyses run after step 06.

**01_Sample_Counts** — Quick local check of sample counts at each filtering stage by parsing QZA archives directly.
```bash
bash code/07_Analysis/01_Diagnostics/01_Sample_Counts/01_diagnostic_sample_counts_local.sh
```
Requires `h5py` (`pip install h5py`). Alternatively, upload `.qzv` files to https://view.qiime2.org.

**02_Assignment_Diagnostics** — Compares taxonomy assignment rates between PS1 and PS2 primer sets, generates diagnostic plots and a markdown report.
```bash
bash code/07_Analysis/01_Diagnostics/02_Assignment_Diagnostics/run_assignment_diagnostics_local.sh
```
Output: `output_data/07_Analysis/assignment_diagnostics/`

**03_Filtering_Loss_Analysis** — Tracks read and feature loss at each pipeline stage to identify where samples drop out.
```bash
bash code/07_Analysis/01_Diagnostics/03_Filtering_Loss_Analysis/run_filtering_loss_analysis_local.sh
```

**04_Local_Barcode_Analysis** — Analyzes the contribution of Schmidt lab's local barcodes to taxonomy assignments.
```bash
bash code/07_Analysis/01_Diagnostics/04_Local_Barcode_Analysis/run_local_barcode_analysis_local.sh
# or on HPCC:
bash code/07_Analysis/01_Diagnostics/04_Local_Barcode_Analysis/run_local_barcode_analysis.sh
```

**05_BLAST_Unassigned** — BLASTs high-abundance unassigned ASVs against NCBI nt to characterize what was not assigned (database gaps, non-target amplification, novel diversity). See `code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/BLAST_Unassigned_README.md` for full details.
```bash
# On HPCC — parallel submission (recommended)
bash code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/blast_unassigned_asvs.sh 100 75
bash code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/blast_unassigned_parallel_submit.sh ps1
bash code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/blast_unassigned_parallel_submit.sh ps2
Rscript code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/analyze_blast_results.R
```

#### 02_Analysis

Ecological analyses and manuscript figures using phyloseq, vegan, and ggplot2. Helper functions are shared across all Rmd files via `gut_analysis_2026_helpers.R`.

Render Rmd files manually in RStudio or with `rmarkdown::render()`. Outputs go to `output_data/07_Analysis/`.

| Script | Description |
|--------|-------------|
| `01_metadata_summary.Rmd` | Sample metadata summary tables (Word document output) |
| `02_bale_creek_comprehensive_analysis.Rmd` | Full community dietary analysis for the Bale Creek mormyrid assemblage: alpha/beta diversity, PERMANOVA, taxonomic composition |
| `03_pk_dietary_variation_complete.Rmd` | Dietary variation in *Paramormyrops kingsleyae*: geographic effects, individual variation, sex differences |
| `04_pk_full_permanova_model.Rmd` | Comprehensive PERMANOVA model for *P. kingsleyae* diet composition |
| `05_pk_simper_selectivity.Rmd` | SIMPER analysis and prey selectivity indices (Ivlev's E, Chesson's Alpha, Jacobs' D) |
| `gut_analysis_2026_helpers.R` | Shared helper functions: data loading, alpha/beta diversity, SIMPER drill-down, selectivity indices, visualization utilities |

## Technical Specifications

### Primer Sets

**Primerset 1 (ps1) - Mliniprobe primers:**
- Forward: `GGWACWGGWTGAACWGTWTAYCCYCC`
- Reverse: `TCDGGRTGNCCRAARAAYCA`
- Truncation: F=224, R=230

**Primerset 2 (ps2) - Leray primers:**
- Forward: `GCHCCHGAYATRGCHTTYCC`
- Reverse: `ARYATDGTRATDGCHCCDGC`
- Truncation: F=230, R=221

### Reference Database Composition

- **COInr (2022-05-06):** Filtered to Metazoa, minimum family-level taxonomy
- **Schmidt barcodes:** 429 vouchered invertebrate sequences from Gabon
- **Region:** Restricted to Leray primer amplicon (280–345 bp)
- **Format:** QIIME2 artifacts for classifier training

### SLURM Job Submission

Standard pattern for submitting jobs to HPCC:

```bash
sbatch --job-name <name> \
       --output ${root}/output_data/slurm_logs/<logfile> \
       --export=root=${root},var1=val1 \
       script.sb
```

- Shell scripts (`.sh`) iterate over runs/primers and submit batch jobs
- SLURM batch scripts (`.sb`) contain resource requirements and singularity commands
- Logs written to `output_data/slurm_logs/`

### Singularity Execution Patterns

**QIIME2 commands:**
```bash
singularity exec --bind $TMPDIR:/home/qiime2/q2cli $qiime_image \
  qiime <plugin> <action> [options]
```

**mkCOInr Perl scripts:**
```bash
singularity exec ${mkcoinr_image} \
  perl ${mkcoinr_path}/<script.pl> [options]
```

**R scripts (local or HPCC):**
```bash
singularity exec ${rimage} Rscript <script.R> [args]
# or locally:
Rscript <script.R> [args]
```

### Resource Requirements

- **DADA2 denoising:** 256GB RAM
- **Classifier training:** 500GB RAM
- **Standard jobs:** 64–128GB RAM
- **Clustering/collapse:** 12GB RAM
- **Analysis scripts (07):** Run locally, no HPCC required

## Data Organization

```
gut_content_2025/
├── gut_contents.env                 # Environment variables (sourced by all scripts)
├── input_data/
│   ├── 01_Filter_And_Denoise/       # Raw FASTQ files
│   ├── 04_Reference_Database/       # COInr, Schmidt barcodes, taxa whitelist
│   └── metadata/                    # Sample metadata
│       ├── wholestudy_metadata.xlsx
│       ├── fish_id_metadata.tsv     # Fish-level metadata
│       ├── primerset1_metadata_final.tsv
│       ├── primerset2_metadata_final.tsv
│       ├── primerset1_sample_to_fish.tsv  # Maps PCR samples to fish IDs
│       └── primerset2_sample_to_fish.tsv
├── output_data/
│   ├── 01_Filter_And_Denoise/       # QIIME2 artifacts from denoising
│   ├── 02_Merge_Data/               # Merged feature tables
│   ├── 03_Clustered_Data/           # Clustered and fish-collapsed ASVs
│   ├── 04_Reference_Database/       # Trained classifiers
│   ├── 05_Taxanomic_Assignment/     # Taxonomy tables (vsearch + naive Bayes)
│   ├── 06_Generate_Output/          # Final filtered, rarefied tables for R
│   ├── 07_Analysis/                 # Diagnostic reports and analysis outputs
│   └── slurm_logs/                  # SLURM job logs
├── code/                            # Pipeline scripts (01–07)
└── images/                          # Container definitions (Dockerfiles, .def files)
```

## Troubleshooting

### Job Failures

Check SLURM logs in `output_data/slurm_logs/` for error messages.

### Container Issues

- Verify container paths in `gut_contents.env` are correct
- Ensure containers built successfully (check `.sif` file exists)
- Remember: Singularity builds must run on HPCC, not locally

### Memory Errors

- DADA2 requires 256GB RAM
- Classifier training requires 500GB RAM
- Adjust SLURM memory requests in `.sb` files if needed

### Rarefaction Depth

- Check `{primer}_all_p985_table_filtd_NO_HOST.qzv` at https://view.qiime2.org after step 06
- Adjust `--p-sampling-depth` in `code/06_Generate_Output/02_generate_output_host_filtered.sh` to a depth that retains most samples with adequate coverage

### mkCOInr Errors

- Scripts expect specific TSV formats
- Check intermediate outputs for formatting issues
- Verify taxonomy strings follow expected hierarchy

### Missing Samples

- Check read count filtering thresholds (`ps1_lo_thresh`, `ps2_lo_thresh`) in `gut_contents.env`
- Review filtering logs in step 02_MergeData
- Verify sample metadata matches SampleIDs in QIIME2 artifacts
- Run `07_Analysis/01_Diagnostics` scripts to trace where samples are lost

### BLAST Issues

- Verify HPCC BLAST database path: `ls /mnt/research/common-data/Bio/blastdb/v5/nt*`
- See `code/07_Analysis/01_Diagnostics/05_BLAST_unassigned/BLAST_Unassigned_README.md` for full troubleshooting

## Environment Setup

All scripts source `gut_contents.env` which defines critical variables:

```bash
source ${root}/gut_contents.env
```

This file contains:
- Container paths
- Primer sequences and truncation lengths
- Read count thresholds
- Run and primer set definitions
