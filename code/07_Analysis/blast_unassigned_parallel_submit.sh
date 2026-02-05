#!/bin/bash
# Main script to split FASTA and submit parallel BLAST job array
#
# Usage: bash blast_unassigned_parallel_submit.sh <PRIMERSET> [sequences_per_job]
#
# Arguments:
#   PRIMERSET: Either "ps1" or "ps2"
#   sequences_per_job: Number of sequences per array job (default: 5)

set -e

# Get repository root (assuming script is in code/07_Analysis/)
root=$(cd "$(dirname "$0")/../.." && pwd)

# Parse arguments
if [ $# -lt 1 ]; then
    echo "Usage: bash blast_unassigned_parallel_submit.sh <PRIMERSET> [sequences_per_job]"
    echo "  PRIMERSET: ps1 or ps2"
    echo "  sequences_per_job: Number of sequences per job (default: 5)"
    exit 1
fi

PRIMERSET=$1
SEQS_PER_JOB=${2:-5}

# Set up directories
BLAST_DIR="${root}/output_data/07_Analysis/blast_unassigned"
SEQ_DIR="${BLAST_DIR}/sequences"
RESULTS_DIR="${BLAST_DIR}/results"
CHUNKS_DIR="${BLAST_DIR}/chunks/${PRIMERSET}"
TMP_RESULTS="${RESULTS_DIR}/tmp_${PRIMERSET}"

mkdir -p ${CHUNKS_DIR}
mkdir -p ${TMP_RESULTS}

# Input file
INPUT_FASTA="${SEQ_DIR}/${PRIMERSET}_unassigned_highAbund.fasta"

echo "========================================"
echo "Parallel BLAST Array Job Setup"
echo "========================================"
echo "Primerset: ${PRIMERSET}"
echo "Input: ${INPUT_FASTA}"
echo "Sequences per job: ${SEQS_PER_JOB}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FASTA}"
    echo "Please run blast_unassigned_asvs.sh first"
    exit 1
fi

# Count sequences
N_SEQS=$(grep -c "^>" ${INPUT_FASTA})
echo "Total sequences: ${N_SEQS}"

# Calculate number of chunks
N_CHUNKS=$(( (N_SEQS + SEQS_PER_JOB - 1) / SEQS_PER_JOB ))
echo "Number of array jobs: ${N_CHUNKS}"
echo ""

# Split FASTA file into chunks
echo "Splitting FASTA file into chunks..."
awk -v seqs_per_job=${SEQS_PER_JOB} -v outdir="${CHUNKS_DIR}" '
BEGIN {
    chunk = 1
    seq_count = 0
    outfile = sprintf("%s/chunk_%04d.fasta", outdir, chunk)
}
/^>/ {
    if (seq_count > 0 && seq_count % seqs_per_job == 0) {
        close(outfile)
        chunk++
        outfile = sprintf("%s/chunk_%04d.fasta", outdir, chunk)
    }
    seq_count++
}
{
    print > outfile
}
END {
    close(outfile)
}
' ${INPUT_FASTA}

echo "Created ${N_CHUNKS} chunk files in ${CHUNKS_DIR}"
echo ""

# Create logs directory if needed
LOGS_DIR="${root}/output_data/slurm_logs"
mkdir -p ${LOGS_DIR}

# Submit SLURM job array
echo "Submitting SLURM job array..."
JOB_ID=$(sbatch --parsable \
    --array=1-${N_CHUNKS} \
    --output=${LOGS_DIR}/blast_${PRIMERSET}_array-%A_%a.out \
    --export=root=${root},PRIMERSET=${PRIMERSET},N_CHUNKS=${N_CHUNKS} \
    ${root}/code/07_Analysis/blast_unassigned_array.sb)

echo "Job array submitted: ${JOB_ID}"
echo ""
echo "Monitor progress with: squeue -u \$USER -j ${JOB_ID}"
echo "Check logs in: ${LOGS_DIR}/blast_${PRIMERSET}_array-${JOB_ID}_*.out"
echo ""
echo "Final results will be in:"
echo "  ${RESULTS_DIR}/${PRIMERSET}_blast_results.txt"
echo "  ${RESULTS_DIR}/${PRIMERSET}_blast_summary.tsv"
