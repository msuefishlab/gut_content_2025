#!/bin/bash
# ==============================================================================
# Local Barcode Contribution Analysis - HPCC Execution Wrapper
# ==============================================================================
# Purpose: Run the local barcode contribution analysis on HPCC using Singularity
#
# Usage: bash code/07_Analysis/run_local_barcode_analysis.sh
# ==============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get script directory and repository root
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export root="$( cd "${script_dir}/../.." && pwd )"

# Source environment variables
source "${root}/gut_contents.env"

# Verify Singularity image exists
if [ ! -f "${rimage}" ]; then
    echo "ERROR: R Singularity image not found at ${rimage}"
    echo "Please build the image first using the appropriate build script."
    exit 1
fi

# Define output directory
output_dir="${root}/output_data/07_Analysis/local_barcode_contribution"

# Create output directories
mkdir -p "${output_dir}/tables"
mkdir -p "${output_dir}/figures"

# Print header
echo "═══════════════════════════════════════════════════════════════════════"
echo "LOCAL BARCODE CONTRIBUTION ANALYSIS (HPCC Execution)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""
echo "Repository root: ${root}"
echo "Singularity image: ${rimage}"
echo "Output directory: ${output_dir}"
echo ""

# Run the analysis using Singularity
echo "Starting analysis with Singularity..."
echo "─────────────────────────────────────────────────────────────────────"
echo ""

singularity exec "${rimage}" Rscript "${root}/code/07_Analysis/local_barcode_contribution_analysis.R"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "═══════════════════════════════════════════════════════════════════════"
    echo "ANALYSIS COMPLETED SUCCESSFULLY"
    echo "═══════════════════════════════════════════════════════════════════════"
    echo ""
    echo "Review the comprehensive report:"
    echo "  ${output_dir}/local_barcode_contribution_report.md"
    echo ""
    echo "Output files:"
    echo "  Tables: ${output_dir}/tables/"
    echo "  Figures: ${output_dir}/figures/"
    echo ""
else
    echo ""
    echo "═══════════════════════════════════════════════════════════════════════"
    echo "ERROR: Analysis failed"
    echo "═══════════════════════════════════════════════════════════════════════"
    echo ""
    exit 1
fi
