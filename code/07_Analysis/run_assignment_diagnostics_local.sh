#!/bin/bash

# ==============================================================================
# Taxonomy Assignment Diagnostics Analysis - Local Wrapper Script
# ==============================================================================
# Purpose: Execute comprehensive analysis of taxonomy assignment rates between
#          PS1 and PS2 primer sets (LOCAL EXECUTION)
# ==============================================================================

# Determine repository root (script is in code/07_Analysis/)
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export root="$( cd "${script_dir}/../.." && pwd )"

# Create output directory
output_dir=${root}/output_data/07_Analysis/assignment_diagnostics
mkdir -p ${output_dir}/{tables,sequences,figures}

echo "========================================"
echo "Taxonomy Assignment Diagnostics Analysis"
echo "========================================"
echo ""
echo "Repository root: ${root}"
echo "Output directory: ${output_dir}"
echo "Running locally with R (no Singularity)"
echo ""

# Run R analysis script directly
echo "Starting R analysis..."
Rscript ${root}/code/07_Analysis/assignment_diagnostics_analysis.R

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
    echo "  - assignment_diagnostics_report.md (full report)"
    echo "  - figures/ (all diagnostic plots)"
    echo "  - sequences/ (FASTA files for BLAST verification)"
    echo "  - tables/ (all data tables)"
    echo ""
else
    echo ""
    echo "ERROR: Analysis failed. Check error messages above."
    exit 1
fi
