#!/bin/bash
# ==============================================================================
# Local Barcode Contribution Analysis - Local Execution Wrapper
# ==============================================================================
# Purpose: Run the local barcode contribution analysis on a local machine
#          (without Singularity)
#
# Usage: bash code/07_Analysis/run_local_barcode_analysis_local.sh
# ==============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get script directory and repository root
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export root="$( cd "${script_dir}/../../../.." && pwd )"

# Define output directory
output_dir="${root}/output_data/07_Analysis/local_barcode_contribution"

# Create output directories
mkdir -p "${output_dir}/tables"
mkdir -p "${output_dir}/figures"

# Print header
echo "═══════════════════════════════════════════════════════════════════════"
echo "LOCAL BARCODE CONTRIBUTION ANALYSIS (Local Execution)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""
echo "Repository root: ${root}"
echo "Output directory: ${output_dir}"
echo ""

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R."
    exit 1
fi

# Check for required R packages
echo "Checking R dependencies..."
Rscript -e "
required_packages <- c('tidyverse', 'qiime2R', 'scales')
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat(sprintf('ERROR: Missing R packages: %s\n', paste(missing, collapse = ', ')))
  cat('Install with: install.packages(c(\"%s\"))\n', paste(missing, collapse = '\", \"'))
  quit(status = 1)
} else {
  cat('✓ All required packages found\n')
}
"

if [ $? -ne 0 ]; then
    exit 1
fi

echo ""

# Run the analysis
echo "Starting analysis..."
echo "─────────────────────────────────────────────────────────────────────"
echo ""

Rscript "${root}/code/07_Analysis/01_Diagnostics/04_Local_Barcode_Analysis/local_barcode_contribution_analysis.R"

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
