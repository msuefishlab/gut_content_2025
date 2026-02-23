#!/bin/bash

# ==============================================================================
# Filtering Loss Analysis - Local Wrapper Script
# ==============================================================================
# Purpose: Quantify taxonomy filtering loss to determine if differential
#          filtering rates could explain compositional differences (LOCAL)
# ==============================================================================

# Determine repository root (script is in code/07_Analysis/01_Diagnostics/03_Filtering_Loss_Analysis/)
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export root="$( cd "${script_dir}/../../../.." && pwd )"

# Create output directory
output_dir="${root}/output_data/07_Analysis/filtering_loss_analysis"
mkdir -p "${output_dir}"/{tables,figures,report}

echo "========================================"
echo "Filtering Loss Analysis"
echo "========================================"
echo ""
echo "Repository root: ${root}"
echo "Output directory: ${output_dir}"
echo "Running locally with R (no Singularity)"
echo ""

# Run R analysis script directly
echo "Starting R analysis..."
Rscript "${root}/code/07_Analysis/01_Diagnostics/03_Filtering_Loss_Analysis/filtering_loss_analysis.R"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Analysis completed successfully!"
    echo "========================================"
    echo ""
    echo "Review outputs at: ${output_dir}"
    echo ""
    echo "Key files:"
    echo "  - report/filtering_loss_report.md (full report)"
    echo "  - figures/ (5 PDF plots)"
    echo "  - tables/ (6 CSV files)"
    echo ""
    echo "Quick summary:"
    echo "  - global_filtering_summary.csv (overall loss rates)"
    echo "  - statistical_tests_summary.csv (significance tests)"
    echo ""
else
    echo ""
    echo "ERROR: Analysis failed. Check error messages above."
    exit 1
fi
