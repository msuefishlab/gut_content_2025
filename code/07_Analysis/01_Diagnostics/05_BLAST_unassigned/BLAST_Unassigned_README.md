# BLAST Analysis of Unassigned ASVs

This directory contains scripts to characterize unassigned ASVs from the gut content metabarcoding analysis by BLASTing them against the NCBI nt database.

## Overview

The pipeline extracts high-abundance unassigned ASVs from each primer set and BLASTs them against NCBI nt to determine:
- Are they metazoan sequences missing from the reference database?
- Are they non-target amplifications (bacteria, fungi, plants)?
- Are they poor-quality sequences?
- Are they novel/divergent lineages?

## Scripts

### 1. `blast_unassigned_asvs.sh`

Extracts and prepares sequences for BLAST analysis.

**Usage:**
```bash
source gut_contents.env
bash code/07_Analysis/blast_unassigned_asvs.sh [min_reads] [n_samples]
```

**Parameters:**
- `min_reads` (default: 100): Minimum read count to consider "high-abundance"
- `n_samples` (default: 75): Number of ASVs to randomly sample per primer set

**Output:**
- `output_data/07_Analysis/blast_unassigned/sequences/`
  - `ps1_unassigned_highAbund.fasta`: FASTA file for PS1
  - `ps2_unassigned_highAbund.fasta`: FASTA file for PS2
  - `ps1_unassigned_metadata.tsv`: Metadata for PS1 ASVs
  - `ps2_unassigned_metadata.tsv`: Metadata for PS2 ASVs

**Notes:**
- Sampling is weighted by abundance to prioritize high-read ASVs
- Sequences are exported from QIIME2 artifacts if not already done
- ASV names include read counts for reference

### 2a. `blast_unassigned_parallel_submit.sh` (RECOMMENDED)

**Parallel version using SLURM job arrays for faster execution.**

**Usage:**
```bash
# On HPCC
source gut_contents.env
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1 [sequences_per_job]
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps2 [sequences_per_job]
```

**Parameters:**
- `PRIMERSET`: Either "ps1" or "ps2"
- `sequences_per_job` (optional, default: 5): Number of sequences per array job

**Resources (per array task):**
- Time: 1 hour
- CPUs: 2
- Memory: 4GB

**Advantages:**
- **10-100x faster wall time** - runs many BLAST jobs in parallel
- **Better resource utilization** - small jobs complete quickly
- **Automatic merging** - results combined automatically
- **Fault tolerance** - individual job failures don't affect others

**Output:**
- Same as serial version (see below)
- Temporary files automatically cleaned up

**How it works:**
1. Splits input FASTA into chunks (5 sequences per chunk by default)
2. Submits SLURM job array with one task per chunk
3. Each task BLASTs its chunk independently
4. Final task automatically merges results and generates summary

### 2b. `blast_unassigned_submit.sb` (Serial version)

**Single-job BLAST submission (slower but simpler).**

**Usage:**
```bash
# On HPCC
sbatch --export=root=/path/to/repo,PRIMERSET=ps1 \
       --output=output_data/slurm_logs/blast_ps1-%j.out \
       code/07_Analysis/blast_unassigned_submit.sb

sbatch --export=root=/path/to/repo,PRIMERSET=ps2 \
       --output=output_data/slurm_logs/blast_ps2-%j.out \
       code/07_Analysis/blast_unassigned_submit.sb
```

**Parameters:**
- `root`: Path to repository root
- `PRIMERSET`: Either "ps1" or "ps2"

**Resources:**
- Time: 12 hours
- CPUs: 8
- Memory: 32GB

**Note:** Use parallel version for faster results. This serial version is provided as a fallback.

**Output (both versions):**
- `output_data/07_Analysis/blast_unassigned/results/`
  - `ps1_blast_results.txt`: Full BLAST output (tabular format)
  - `ps1_blast_summary.tsv`: Top hit per ASV with taxonomy
  - `ps2_blast_results.txt`: Full BLAST output for PS2
  - `ps2_blast_summary.tsv`: Top hit per ASV for PS2

### 3. `blast_unassigned_run_all.sh`

Wrapper script to run the entire pipeline (prepare + submit).

**Usage:**
```bash
# On HPCC
source gut_contents.env
bash code/07_Analysis/blast_unassigned_run_all.sh [min_reads] [n_samples]
```

This script:
1. Runs sequence preparation (`blast_unassigned_asvs.sh`)
2. Submits BLAST jobs for both primer sets
3. Provides job IDs for monitoring

### 4. `analyze_blast_results.R`

R script to analyze BLAST results and generate summary statistics and visualizations.

**Usage:**
```bash
# After BLAST jobs complete
Rscript code/07_Analysis/analyze_blast_results.R
```

**Output:**
- `output_data/07_Analysis/blast_unassigned/results/`
  - `blast_analysis_summary.csv`: Detailed results with categorization
- `output_data/07_Analysis/blast_unassigned/plots/`
  - `hit_quality.pdf`: Hit quality distribution
  - `taxonomic_composition.pdf`: Taxonomic breakdown
  - `identity_distribution.pdf`: BLAST identity histograms
  - `abundance_vs_identity.pdf`: Scatter plot of reads vs identity
  - `blast_analysis_combined.pdf`: Combined multi-panel figure
- `output_data/07_Analysis/blast_unassigned/`
  - `BLAST_Analysis_Report.md`: Markdown summary report

## Workflow

### Quick Start - Parallel (RECOMMENDED)

```bash
# 1. Set up environment
source gut_contents.env

# 2. Extract and prepare sequences
bash code/07_Analysis/blast_unassigned_asvs.sh 100 75

# 3. Submit parallel BLAST jobs (on HPCC)
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps2

# 4. Monitor jobs (will show array of jobs)
squeue -u $USER

# 5. Once jobs complete, analyze results
Rscript code/07_Analysis/analyze_blast_results.R
```

**Expected time:** Minutes to 1-2 hours (vs 12 hours for serial version)

### Quick Start - Serial (Fallback)

```bash
# 1. Set up environment
source gut_contents.env

# 2. Run full pipeline (on HPCC)
bash code/07_Analysis/blast_unassigned_run_all.sh

# 3. Monitor jobs
squeue -u $USER

# 4. Once jobs complete, analyze results
Rscript code/07_Analysis/analyze_blast_results.R
```

### Step-by-Step (Parallel)

```bash
# 1. Extract and prepare sequences
source gut_contents.env
bash code/07_Analysis/blast_unassigned_asvs.sh 100 75

# 2. Submit parallel BLAST jobs (on HPCC)
# Default: 5 sequences per job
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps2

# Or customize chunk size (e.g., 10 sequences per job)
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1 10

# 3. Check job status (array jobs show as jobID_[1-N])
squeue -u $USER

# 4. Monitor individual array task
tail -f blast_array-JOBID_TASKID.out

# 5. Results are automatically merged when all tasks complete
# 6. Analyze results
Rscript code/07_Analysis/analyze_blast_results.R
```

### Step-by-Step (Serial)

```bash
# 1. Extract and prepare sequences
source gut_contents.env
bash code/07_Analysis/blast_unassigned_asvs.sh 100 75

# 2. Submit BLAST jobs (on HPCC)
sbatch --export=root=${root},PRIMERSET=ps1 \
       --output=output_data/slurm_logs/blast_ps1-%j.out \
       code/07_Analysis/blast_unassigned_submit.sb

sbatch --export=root=${root},PRIMERSET=ps2 \
       --output=output_data/slurm_logs/blast_ps2-%j.out \
       code/07_Analysis/blast_unassigned_submit.sb

# 3. Check job status
squeue -u $USER

# 4. View log files
tail -f output_data/slurm_logs/blast_ps1-*.out

# 5. Analyze results
Rscript code/07_Analysis/analyze_blast_results.R
```

## Customization

### Parallel vs Serial: Which to Use?

**Use Parallel (recommended):**
- Default for most cases - faster and more efficient
- When you have >10 sequences to BLAST
- When you want results quickly

**Use Serial:**
- Very few sequences (<10)
- Testing or debugging
- When HPCC job queue is very busy

### Tune Parallel Performance

The `sequences_per_job` parameter controls parallelization:

```bash
# More parallelization (faster, more jobs)
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1 1

# Balanced (default, recommended for 50-100 sequences)
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1 5

# Less parallelization (fewer jobs, longer per job)
bash code/07_Analysis/blast_unassigned_parallel_submit.sh ps1 10
```

**Guidelines:**
- **1-2 seqs/job**: Maximum parallelization, best for >100 sequences
- **5 seqs/job** (default): Good balance for 50-100 sequences
- **10-20 seqs/job**: Fewer jobs, better for <50 sequences

### Adjust abundance threshold

To focus on even higher abundance ASVs (e.g., ≥500 reads):

```bash
bash code/07_Analysis/blast_unassigned_asvs.sh 500 50
```

### Change sample size

To BLAST more sequences per primer set:

```bash
bash code/07_Analysis/blast_unassigned_asvs.sh 100 150
```

### Modify BLAST parameters

Edit `blast_unassigned_array.sb` (parallel) or `blast_unassigned_submit.sb` (serial) to adjust:
- `-max_target_seqs`: Number of hits per query (default: 10)
- `-evalue`: E-value threshold (default: 1e-5)
- `-num_threads`: CPU threads (default: 2 for parallel, 8 for serial)

## Output Interpretation

### BLAST Summary Columns

- `ASV`: ASV identifier (MD5 hash)
- `Reads`: Total read count for this ASV
- `Top_Hit_Species`: Scientific name from top BLAST hit
- `Top_Hit_Common`: Common name (if available)
- `Percent_Identity`: Sequence identity (%)
- `E_value`: BLAST E-value
- `Bit_Score`: BLAST bit score
- `Alignment_Length`: Length of alignment
- `Top_Hit_Title`: Full sequence description

### Key Questions to Address

1. **Are unassigned ASVs metazoan?**
   - Check `Likely_Metazoan` column in analysis output
   - High proportion suggests reference database gaps

2. **Are they COI sequences?**
   - Check `Likely_COI` column
   - High proportion suggests taxonomy assignment issue, not marker gene issue

3. **What's the identity distribution?**
   - ≥97%: Very close matches, likely same species
   - 90-97%: Good matches, likely congeners
   - 80-90%: Moderate matches, family-level
   - <80%: Poor matches, potential novel diversity

4. **Are they contaminants?**
   - Check for bacterial, fungal, or plant hits
   - May indicate non-target amplification

## Troubleshooting

### No sequences found

If you get "No ASVs found with >= N reads":
- Lower the `min_reads` threshold
- Check that `removed_asvs_full_list.csv` contains data

### BLAST database not found

Ensure the HPCC BLAST database path is correct:
```bash
ls /mnt/research/common-data/Bio/blastdb/v5/nt*
```

If path changed, update `blast_unassigned_submit.sb`.

### Out of memory

If BLAST runs out of memory:
- Reduce `n_samples` to BLAST fewer sequences
- Increase `--mem` in the SBATCH directives

### Analysis script fails

Ensure required R packages are installed:
```R
install.packages(c("tidyverse", "patchwork", "Biostrings"))
```

## Notes

- **Parallel version is 10-100x faster** than serial (minutes to hours vs 12 hours)
- Parallel version automatically cleans up temporary chunk files after merging
- BLAST results are identical between serial and parallel versions
- Results are saved in tabular format for easy parsing
- Random sampling uses set.seed(42) for reproducibility
- Weighted sampling prioritizes high-abundance ASVs
- All scripts follow the project convention of sourcing `gut_contents.env`
- Job array logs are saved as `blast_array-JOBID_TASKID.out` in submission directory

## Citation

If using these scripts, cite the BLAST+ tool:
- Camacho C. et al. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421.
