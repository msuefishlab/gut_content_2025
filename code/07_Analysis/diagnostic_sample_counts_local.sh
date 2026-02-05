#!/bin/bash
#
# Quick local diagnostic (no Singularity needed)
# Counts samples by examining the QZA archive structure
#

root="$(git rev-parse --show-toplevel)"

echo "=========================================="
echo "Quick Sample Count (Local - No Singularity)"
echo "=========================================="
echo ""

for primer in primerset1 primerset2; do
    echo "==================== $primer ===================="
    echo ""

    # Stage 3: After filtering features
    echo "STAGE 3: After taxonomy filter"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd.qza" ]; then
        unzip -p "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd.qza" \
            '*/data/feature-table.biom' 2>/dev/null | \
            python3 -c "import sys, json; d=json.load(sys.stdin); print(f'  Samples: {len(d[\"columns\"])}')" 2>/dev/null || echo "  (Install python3 to parse)"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 4: After host filtering
    echo "STAGE 4: After host filtering (NO_HOST)"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" ]; then
        unzip -p "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" \
            '*/data/feature-table.biom' 2>/dev/null | \
            python3 -c "
import sys, json
d = json.load(sys.stdin)
samples = d['columns']
data_matrix = d['data']

# Calculate per-sample read counts
sample_counts = {}
for i, sample in enumerate(samples):
    sample_id = sample['id']
    # Sum reads across all features for this sample
    total = sum(row[i] for row in data_matrix)
    sample_counts[sample_id] = total

counts = list(sample_counts.values())
counts.sort()

print(f'  Samples: {len(samples)}')
if counts:
    print(f'  Min reads: {min(counts)}')
    print(f'  Median reads: {counts[len(counts)//2]}')
    print(f'  Max reads: {max(counts)}')
    below_5000 = sum(1 for c in counts if c < 5000)
    print(f'  Samples < 5000 reads: {below_5000} (WILL BE DROPPED!)')
" 2>/dev/null || echo "  (Install python3 to parse)"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 5: After rarefaction
    echo "STAGE 5: After rarefaction to 5000"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" ]; then
        unzip -p "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" \
            '*/data/feature-table.biom' 2>/dev/null | \
            python3 -c "import sys, json; d=json.load(sys.stdin); print(f'  Samples: {len(d[\"columns\"])}')" 2>/dev/null || echo "  (Install python3 to parse)"
    else
        echo "  File not found"
    fi
    echo ""
done

echo "=========================================="
echo "RECOMMENDATION:"
echo "  Upload these .qzv files to https://view.qiime2.org"
echo "  to see detailed sample statistics:"
echo ""
echo "  - primerset1_all_p985_table_filtd_NO_HOST.qzv"
echo "  - primerset2_all_p985_table_filtd_NO_HOST.qzv"
echo ""
echo "  Then adjust the rarefaction depth in:"
echo "  code/06_Generate_Output/02_generate_output_host_filtered.sh:72"
echo "=========================================="
