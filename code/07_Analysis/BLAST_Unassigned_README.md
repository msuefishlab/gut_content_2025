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

### 2. `blast_unassigned_submit.sb`

SLURM script to run BLAST on HPCC.

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

**Output:**
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

### Quick Start (Recommended)

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

### Step-by-Step

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

### Adjust abundance threshold

To focus on even higher abundance ASVs (e.g., ≥500 reads):

```bash
bash code/07_Analysis/blast_unassigned_run_all.sh 500 50
```

### Change sample size

To BLAST more sequences per primer set:

```bash
bash code/07_Analysis/blast_unassigned_run_all.sh 100 150
```

### Modify BLAST parameters

Edit `blast_unassigned_submit.sb` to adjust:
- `-max_target_seqs`: Number of hits per query (default: 10)
- `-evalue`: E-value threshold (default: 1e-5)
- `-num_threads`: CPU threads (default: uses all allocated)

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

- BLAST against nt can be slow (several hours per primer set)
- Results are saved in tabular format for easy parsing
- Random sampling uses set.seed(42) for reproducibility
- Weighted sampling prioritizes high-abundance ASVs
- All scripts follow the project convention of sourcing `gut_contents.env`

## Citation

If using these scripts, cite the BLAST+ tool:
- Camacho C. et al. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421.
