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
3. **Cluster Reads** - Cluster ASVs to reduce redundancy
4. **Reference Database** - Build custom COInr + Schmidt database
5. **Taxonomic Assignment** - Assign taxonomy with vsearch/Bayes
6. **Generate Output** - Export final tables and visualizations

## Container Management

**IMPORTANT:** Singularity build commands must be run on HPCC, not locally.

### QIIME2 Container

```bash
cd images/
singularity build --remote qiime2.sif qiime2_mod.def
```

Alternative (QIIME2 2019.1):
```bash
singularity build qiime2019.sif docker://quay.io/qiime2/core:2019.1
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

### QIIME1 Container (legacy)

```bash
singularity build qiime1.sif docker://mbari/qiime1:latest
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
```

All scripts submit SLURM jobs via `.sb` files. Check logs in `output_data/slurm_logs/`.

### 02_MergeData

Filters and merges denoised data across sequencing runs and primer sets.

**What it does:**
- Filters out samples below read count thresholds (defined in `gut_contents.env`)
- Merges feature tables across runs
- Combines primerset1 and primerset2 data

**How to run:**
```bash
bash code/02_MergeData/01_remove_low_reads.sh
```

### 03_Cluster_Reads

Clusters similar ASVs to reduce dataset complexity while maintaining biological signal.

**What it does:**
- Uses QIIME2 vsearch clustering
- Groups similar sequences at defined similarity threshold

**How to run:**
```bash
bash code/03_Cluster_Reads/01_run_cluster.sh
```

### 04_Reference_Database

Creates custom reference database for taxonomic assignment.

**Key workflow:**
1. Download COInr database from Zenodo (2022-05-06 version)
2. Filter COInr to Metazoa at family level minimum
3. Process Schmidt lab's 429 sequenced barcodes from Gabon invertebrates
4. Merge and dereplicate sequences
5. Restrict to Leray primer amplicon region (280-345 bp)
6. Train QIIME2 naive Bayes classifier

**How to run:**
```bash
bash code/04_Reference_Database/01_run_classifier_training.sh
```

See `code/04_Reference_Database/README.md` for detailed workflow documentation.

Uses mkCOInr Perl scripts extensively via Singularity container.

### 05_Taxanomic_Assignment

Assigns taxonomy to clustered ASVs using two complementary approaches.

**What it does:**
- **vsearch:** Exact sequence matching for high-confidence assignments
- **Bayes classifier:** Probabilistic assignment for novel sequences

**How to run:**
```bash
bash code/05_Taxanomic_Assignment/01_run_vsearch.sh
```

### 06_Generate_Output

Exports final tables and creates visualizations for downstream analysis in R.

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
- **Region:** Restricted to Leray primer amplicon (280-345 bp)
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

### Resource Requirements

- **DADA2 denoising:** 256GB RAM
- **Classifier training:** 500GB RAM
- **Standard jobs:** 64-128GB RAM

## Data Organization

```
gut_content_2025/
├── input_data/
│   ├── 01_Filter_And_Denoise/       # Raw FASTQ files
│   ├── 04_Reference_Database/        # COInr, Schmidt barcodes
│   └── metadata/                     # Sample metadata (wholestudy_metadata.xlsx)
├── output_data/
│   ├── 01_Filter_And_Denoise/       # QIIME2 artifacts from denoising
│   ├── 02_MergeData/                # Merged feature tables
│   ├── 03_Cluster_Reads/            # Clustered ASVs
│   ├── 04_Reference_Database/        # Trained classifiers
│   ├── 05_Taxanomic_Assignment/     # Taxonomy tables
│   ├── 06_Generate_Output/          # Final tables for R
│   └── slurm_logs/                  # SLURM job logs
├── code/                            # Pipeline scripts (01-06)
└── images/                          # Singularity containers
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

### mkCOInr Errors

- Scripts expect specific TSV formats
- Check intermediate outputs for formatting issues
- Verify taxonomy strings follow expected hierarchy

### Missing Samples

- Check read count filtering thresholds in `gut_contents.env`
- Review filtering logs in step 02_MergeData
- Verify sample metadata matches SampleIDs in QIIME2 artifacts

## Environment Setup

All scripts source `gut_contents.env` which defines critical variables:

```bash
source ${root}/gut_contents.env
```

This file contains:
- Container paths
- Primer sequences
- Truncation lengths
- Read count thresholds
- Run and primer set definitions
