root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"


root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	metadata=${root}/input_data/metadata/${primer}_metadata_final.tsv
	output_dir=${root}/output_data/06_Generate_Output

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table summarize \
	--i-table ${output_dir}/${primer}_all_p985_table_filtd.qza \
	--o-visualization ${output_dir}/${primer}_all_p985_table_filtd.qzv \
	--m-sample-metadata-file ${metadata}
  
  	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table tabulate-seqs \
  	--i-data ${output_dir}/${primer}_all_p985_seqs_Filtd.qza \
  	--o-visualization ${output_dir}/${primer}_all_p985_seqs_Filtd.qzv

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime taxa barplot \
	--i-table ${output_dir}/${primer}_all_p985_table_filtd.qza --i-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza \
	--m-metadata-file ${metadata} \
	--o-visualization ${output_dir}/${primer}_taxa_bar_plots_no_rarify.qzv

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime diversity alpha-rarefaction \
  	--i-table ${output_dir}/${primer}_all_p985_table_filtd.qza \
  	--p-min-depth 1000 --p-max-depth 12000 \
  	--p-metrics shannon observed_features \
  	--o-visualization ${output_dir}/${primer}_all_p985_rarefactionCurve.qzv

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table rarefy \
  	--i-table  ${output_dir}/${primer}_all_p985_table_filtd.qza \
  	--p-sampling-depth 5000 \
  	--o-rarefied-table  ${output_dir}/${primer}_all_p985_table_filtd.rarefied.qza

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime taxa barplot \
	--i-table ${output_dir}/${primer}_all_p985_table_filtd.rarefied.qza --i-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza \
	--m-metadata-file ${metadata} \
	--o-visualization ${output_dir}/${primer}_all_p985_table_filtd.rarefied.qzv

done

