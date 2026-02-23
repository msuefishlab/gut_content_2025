#!/bin/bash
#
# Quick local diagnostic (no Singularity needed)
# Counts samples by examining the QZA archive structure
#
# NOTE: BIOM files in QIIME2 are in HDF5 format, not JSON.
# This script requires h5py: pip install h5py
#

root="$(git rev-parse --show-toplevel)"

# Check if h5py is available
if ! python3 -c "import h5py" 2>/dev/null; then
    echo "=========================================="
    echo "ERROR: h5py not installed"
    echo "=========================================="
    echo ""
    echo "BIOM files are in HDF5 format and require h5py to parse."
    echo ""
    echo "To install: pip install h5py"
    echo ""
    echo "Or just view the QZV files at https://view.qiime2.org:"
    echo "  - primerset1_all_p985_table_filtd_NO_HOST.qzv"
    echo "  - primerset2_all_p985_table_filtd_NO_HOST.qzv"
    echo "=========================================="
    exit 1
fi

echo "=========================================="
echo "Quick Sample Count (Local - No Singularity)"
echo "=========================================="
echo ""

for primer in primerset1 primerset2; do
    echo "==================== $primer ===================="
    echo ""

    # Stage 3: After filtering features
    echo "STAGE 3: After taxonomy filter (grouped)"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_grouped.qza" ]; then
        # Extract and analyze with h5py
        tmpdir=$(mktemp -d)
        unzip -q "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_grouped.qza" -d "$tmpdir"
        biom_file=$(find "$tmpdir" -name "feature-table.biom")
        python3 <<EOF
import h5py
try:
    with h5py.File('$biom_file', 'r') as f:
        n_samples = f['sample/ids'].shape[0]
        print(f'  Samples: {n_samples}')
except Exception as e:
    print(f'  (Error: {e})')
EOF
        rm -rf "$tmpdir"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 4: After host filtering
    echo "STAGE 4: After host filtering (NO_HOST)"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" ]; then
        tmpdir=$(mktemp -d)
        unzip -q "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.qza" -d "$tmpdir"
        biom_file=$(find "$tmpdir" -name "feature-table.biom")
        python3 <<EOF
import h5py
import numpy as np

try:
    with h5py.File('$biom_file', 'r') as f:
        n_samples = f['sample/ids'].shape[0]

        # Read the sparse matrix data
        indptr = f['sample/matrix/indptr'][:]
        data = f['sample/matrix/data'][:]

        # Calculate per-sample read counts
        sample_counts = []
        for i in range(n_samples):
            start = indptr[i]
            end = indptr[i + 1]
            total = data[start:end].sum()
            sample_counts.append(int(total))

        sample_counts.sort()

        print(f'  Samples: {n_samples}')
        if sample_counts:
            print(f'  Min reads: {min(sample_counts)}')
            print(f'  Median reads: {sample_counts[len(sample_counts)//2]}')
            print(f'  Max reads: {max(sample_counts)}')
            below_5000 = sum(1 for c in sample_counts if c < 5000)
            print(f'  Samples < 5000 reads: {below_5000} (WILL BE DROPPED in rarefaction!)')
except Exception as e:
    print(f'  (Error: {e})')
EOF
        rm -rf "$tmpdir"
    else
        echo "  File not found"
    fi
    echo ""

    # Stage 5: After rarefaction
    echo "STAGE 5: After rarefaction to 5000"
    if [ -f "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" ]; then
        tmpdir=$(mktemp -d)
        unzip -q "${root}/output_data/06_Generate_Output/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza" -d "$tmpdir"
        biom_file=$(find "$tmpdir" -name "feature-table.biom")
        python3 <<EOF
import h5py
try:
    with h5py.File('$biom_file', 'r') as f:
        n_samples = f['sample/ids'].shape[0]
        print(f'  Samples: {n_samples}')
except Exception as e:
    print(f'  (Error: {e})')
EOF
        rm -rf "$tmpdir"
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
