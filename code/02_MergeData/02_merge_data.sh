#!/bin/bash --login

root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

input_data=${root}/output_data/01_Filter_And_Denoise
output_data=${root}/output_data/02_Merge_Data

mkdir -p ${output_data}

mkdir -p ${SCRATCH}/tmp
export TMPDIR=${SCRATCH}/tmp

#primerset1
singularity exec --bind $TMPDIR:/home/qiime2/q2cli $qiime_image qiime feature-table merge \
  --i-tables ${output_data}/main_run_primerset1_filtered_table.qza \
  --i-tables ${input_data}/repeat_run_primerset1_denoised/table.qza  \
  --o-merged-table ${output_data}/primerset1_all_table.qza


singularity exec --bind $TMPDIR:/home/qiime2/q2cli $qiime_image qiime feature-table merge-seqs \
  --i-data ${output_data}/main_run_primerset1_filtered_rep_seqs.qza \
  --i-data ${input_data}/repeat_run_primerset1_denoised/representative_sequences.qza  \
  --o-merged-data ${output_data}/primerset1_all_rep_seqs.qza

#primerset2
singularity exec --bind $TMPDIR:/home/qiime2/q2cli $qiime_image qiime feature-table merge \
  --i-tables ${output_data}/main_run_primerset2_filtered_table.qza \
  --i-tables ${input_data}/repeat_run_primerset2_denoised/table.qza  \
  --o-merged-table ${output_data}/primerset2_all_table.qza


singularity exec --bind $TMPDIR:/home/qiime2/q2cli $qiime_image qiime feature-table merge-seqs \
  --i-data ${output_data}/main_run_primerset2_filtered_rep_seqs.qza \
  --i-data ${input_data}/repeat_run_primerset2_denoised/representative_sequences.qza  \
  --o-merged-data ${output_data}/primerset2_all_rep_seqs.qza