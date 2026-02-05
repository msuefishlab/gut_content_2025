#!/bin/bash
# run_filtering_loss_analysis.sh
#
# Execute filtering loss analysis using R singularity container
#
# Usage: bash code/07_Analysis/run_filtering_loss_analysis.sh

# Source environment variables
source "${root}/gut_contents.env"

# Run analysis
echo "========================================="
echo "Filtering Loss Analysis"
echo "========================================="
echo ""
echo "Running R script in singularity container..."
echo ""

singularity exec "${rimage}" Rscript "${root}/code/07_Analysis/filtering_loss_analysis.R"

echo ""
echo "========================================="
echo "Analysis complete!"
echo "Results saved to: output_data/07_Analysis/filtering_loss_analysis/"
echo "========================================="
