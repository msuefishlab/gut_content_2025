#!/bin/bash
#
# Modified pipeline: Filter host (Chordata) sequences BEFORE rarefaction
# This corrects for differential host DNA amplification between waveform types
#

root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"
mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
    echo "========================================"
    echo "Processing: $primer"
    echo "========================================"
    
    metadata=${root}/input_data/metadata/${primer}_metadata_final.tsv
    output_dir=${root}/output_data/06_Generate_Output
    
    # --- STEP 1: Summarize UNFILTERED table (for comparison) ---
    echo "Step 1: Summarizing unfiltered table..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table summarize \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd.qza \
        --o-visualization ${output_dir}/${primer}_all_p985_table_filtd.qzv \
        --m-sample-metadata-file ${metadata}
    
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table tabulate-seqs \
        --i-data ${output_dir}/${primer}_all_p985_seqs_Filtd.qza \
        --o-visualization ${output_dir}/${primer}_all_p985_seqs_Filtd.qzv
    
    # --- STEP 2: Filter out host (Chordata) sequences ---
    echo "Step 2: Filtering host (Chordata) sequences..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime taxa filter-table \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd.qza \
        --i-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza \
        --p-exclude "Chordata,Vertebrata,Actinopterygii,Craniata" \
        --p-mode "contains" \
        --o-filtered-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qza
    
    # --- STEP 3: Summarize HOST-FILTERED table ---
    echo "Step 3: Summarizing host-filtered table..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table summarize \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qza \
        --o-visualization ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qzv \
        --m-sample-metadata-file ${metadata}
    
    # --- STEP 4: Taxa barplot (prey only, unrarefied) ---
    echo "Step 4: Generating taxa barplot (prey only, unrarefied)..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime taxa barplot \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qza \
        --i-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza \
        --m-metadata-file ${metadata} \
        --o-visualization ${output_dir}/${primer}_taxa_bar_plots_NO_HOST_no_rarefy.qzv
    
    # --- STEP 5: Rarefaction curves (prey only) ---
    # IMPORTANT: Check the _NO_HOST.qzv summary to determine appropriate depth range
    # You may need to adjust --p-min-depth and --p-max-depth based on prey-only read counts
    echo "Step 5: Generating rarefaction curves (prey only)..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime diversity alpha-rarefaction \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qza \
        --p-min-depth 500 \
        --p-max-depth 8000 \
        --p-metrics shannon observed_features \
        --o-visualization ${output_dir}/${primer}_all_p985_rarefactionCurve_NO_HOST.qzv
    
    # --- STEP 6: Rarefy prey-only table ---
    # IMPORTANT: Adjust --p-sampling-depth based on the _NO_HOST.qzv summary!
    # Choose a depth that retains most samples while providing adequate coverage
    # This may need to be LOWER than 5000 since host reads are now removed
    echo "Step 6: Rarefying prey-only table..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table rarefy \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.qza \
        --p-sampling-depth 3000 \
        --o-rarefied-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza
    
    # --- STEP 7: Taxa barplot (prey only, rarefied) ---
    echo "Step 7: Generating taxa barplot (prey only, rarefied)..."
    singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime taxa barplot \
        --i-table ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qza \
        --i-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza \
        --m-metadata-file ${metadata} \
        --o-visualization ${output_dir}/${primer}_all_p985_table_filtd_NO_HOST.rarefied.qzv
    
    echo "Completed: $primer"
    echo ""
done

echo "========================================"
echo "IMPORTANT NEXT STEPS:"
echo "========================================"
echo ""
echo "1. Check the *_NO_HOST.qzv files to see prey-only read distributions"
echo "2. Adjust rarefaction depth (currently set to 3000) if needed"
echo "   - For PS1: ~70% of reads were prey in BP, ~46% in TP"
echo "   - Original 5000 reads → expect ~2300-3500 prey reads"
echo "3. Re-run cross-validation with the new *_NO_HOST.rarefied.qza files"
echo ""
