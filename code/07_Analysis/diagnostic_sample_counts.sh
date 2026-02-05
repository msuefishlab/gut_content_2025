#!/bin/bash
#
# Diagnostic script to track sample counts through the pipeline
# This helps identify where samples are being lost
#

root="$(git rev-parse --show-toplevel)"
source "${root}/gut_contents.env"

# Create temp directory for exports
temp_dir="${root}/output_data/temp_diagnostic"
mkdir -p "${temp_dir}"

echo "=========================================="
echo "Sample Count Tracking Through Pipeline"
echo "=========================================="
echo ""

for primer in "${primers[@]}"; do
    echo "==================== $primer ===================="
    echo ""

    # Stage 1: After merging (02_MergeData)
    echo "STAGE 1: After merging runs"
    if [ -f "${root}/output_data/02_Merge_Data/${primer}_all_table.qza" ]; then
        singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image \
            qiime tools export \
            --input-path "${root}/output_data/02_Merge_Data/${primer}_all_table.qza" \
            --output-path "${temp_dir}/${primer}_stage1"

        biom summarize-table -i "${temp_dir}/${primer}_stage1/feature-table.biom" > "${temp_dir}/${primer}_stage1_summary.txt"
        samples=$(grep "Num samples:" "${temp_dir}/${primer}_stage1_summary.txt" | awk '{print $3}')
        echo "  Samples: $samples"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 2: After clustering (03_Clustered_Data)
    echo "STAGE 2: After clustering at 98.5%"
    if [ -f "${root}/output_data/03_Clustered_Data/${primer}_all_p985_table.qza" ]; then
        singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image \
            qiime tools export \
            --input-path "${root}/output_data/03_Clustered_Data/${primer}_all_p985_table.qza" \
            --output-path "${temp_dir}/${primer}_stage2"

        biom summarize-table -i "${temp_dir}/${primer}_stage2/feature-table.biom" > "${temp_dir}/${primer}_stage2_summary.txt"
        samples=$(grep "Num samples:" "${temp_dir}/${primer}_stage2_summary.txt" | awk '{print $3}')
        echo "  Samples: $samples"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 3: After filtering features (06_Generate_Output)
    echo "STAGE 3: After filtering features (taxonomy filter)"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd.qza" ]; then
        singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image \
            qiime tools export \
            --input-path "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd.qza" \
            --output-path "${temp_dir}/${primer}_stage3"

        biom summarize-table -i "${temp_dir}/${primer}_stage3/feature-table.biom" > "${temp_dir}/${primer}_stage3_summary.txt"
        samples=$(grep "Num samples:" "${temp_dir}/${primer}_stage3_summary.txt" | awk '{print $3}')
        echo "  Samples: $samples"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 4: After host filtering
    echo "STAGE 4: After removing host sequences"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" ]; then
        singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image \
            qiime tools export \
            --input-path "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" \
            --output-path "${temp_dir}/${primer}_stage4"

        biom summarize-table -i "${temp_dir}/${primer}_stage4/feature-table.biom" > "${temp_dir}/${primer}_stage4_summary.txt"
        samples=$(grep "Num samples:" "${temp_dir}/${primer}_stage4_summary.txt" | awk '{print $3}')
        min_reads=$(grep "Min:" "${temp_dir}/${primer}_stage4_summary.txt" | awk '{print $2}' | sed 's/\..*$//')
        max_reads=$(grep "Max:" "${temp_dir}/${primer}_stage4_summary.txt" | awk '{print $2}' | sed 's/\..*$//')
        median_reads=$(grep "Median:" "${temp_dir}/${primer}_stage4_summary.txt" | awk '{print $2}' | sed 's/\..*$//')

        # Count samples with < 5000 reads
        samples_below_5000=$(biom summarize-table -i "${temp_dir}/${primer}_stage4/feature-table.biom" --observations | \
            tail -n +17 | awk '{print $NF}' | awk '{if ($1 < 5000) count++} END {print count+0}')

        echo "  Samples: $samples"
        echo "  Min reads per sample: $min_reads"
        echo "  Median reads per sample: $median_reads"
        echo "  Max reads per sample: $max_reads"
        echo "  Samples with < 5000 reads (will be dropped by rarefaction): $samples_below_5000"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 5: After rarefaction
    echo "STAGE 5: After rarefaction to 5000 reads"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" ]; then
        singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image \
            qiime tools export \
            --input-path "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" \
            --output-path "${temp_dir}/${primer}_stage5"

        biom summarize-table -i "${temp_dir}/${primer}_stage5/feature-table.biom" > "${temp_dir}/${primer}_stage5_summary.txt"
        samples=$(grep "Num samples:" "${temp_dir}/${primer}_stage5_summary.txt" | awk '{print $3}')
        echo "  Samples: $samples"
    else
        echo "  File not found"
    fi
    echo ""
    echo ""
done

echo "=========================================="
echo "Summary saved to: ${temp_dir}"
echo ""
echo "KEY INSIGHT:"
echo "  Samples are dropped during rarefaction if they have"
echo "  fewer reads than the rarefaction depth (currently 5000)."
echo ""
echo "SOLUTIONS:"
echo "  1. Lower the rarefaction depth in"
echo "     code/06_Generate_Output/02_generate_output_host_filtered.sh:72"
echo "  2. Use the rarefaction curve visualizations to choose"
echo "     an appropriate depth that balances sample retention"
echo "     vs. sequencing depth"
echo "=========================================="
