# Filtering Loss Analysis - Implementation Summary

## Date: 2026-02-04

## Files Created

1. **`filtering_loss_analysis.R`** - Main analysis script (900+ lines)
2. **`run_filtering_loss_analysis.sh`** - Execution wrapper for HPCC
3. **`README_filtering_loss.md`** - Complete documentation

## Next Steps

Run on HPCC:
```bash
source gut_contents.env
bash code/07_Analysis/run_filtering_loss_analysis.sh
```

Results will be saved to: `output_data/07_Analysis/filtering_loss_analysis/`

## What This Analyzes

Quantifies ASVs and reads lost during Family-level taxonomic filtering to determine if differential filtering between PS1 and PS2 could explain compositional differences.

Expected runtime: 2-3 minutes
