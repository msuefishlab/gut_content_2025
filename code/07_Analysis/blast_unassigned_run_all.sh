#!/bin/bash
#
# Wrapper script to prepare and submit BLAST jobs for both primer sets
#
# Usage: bash blast_unassigned_run_all.sh [min_reads] [n_samples]

source ${root}/gut_contents.env

MIN_READS=${1:-100}
N_SAMPLES=${2:-75}

echo "=================================================="
echo "BLAST Analysis of Unassigned ASVs - Full Pipeline"
echo "=================================================="
echo ""

# Step 1: Extract and prepare sequences
echo "Step 1: Extracting and preparing sequences..."
bash ${root}/code/07_Analysis/blast_unassigned_asvs.sh ${MIN_READS} ${N_SAMPLES}

if [ $? -ne 0 ]; then
    echo "ERROR: Sequence preparation failed"
    exit 1
fi

echo ""
echo "Step 2: Submitting BLAST jobs to HPCC..."

# Submit BLAST jobs for both primer sets
JOB1=$(sbatch --parsable \
    --export=root=${root},PRIMERSET=ps1 \
    --output=${root}/output_data/slurm_logs/blast_ps1-%j.out \
    ${root}/code/07_Analysis/blast_unassigned_submit.sb)

JOB2=$(sbatch --parsable \
    --export=root=${root},PRIMERSET=ps2 \
    --output=${root}/output_data/slurm_logs/blast_ps2-%j.out \
    ${root}/code/07_Analysis/blast_unassigned_submit.sb)

echo "  Submitted PS1 BLAST job: ${JOB1}"
echo "  Submitted PS2 BLAST job: ${JOB2}"
echo ""
echo "Monitor jobs with: squeue -u ${USER}"
echo ""
echo "Once complete, run analysis with:"
echo "  Rscript ${root}/code/07_Analysis/analyze_blast_results.R"
echo ""
