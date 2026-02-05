root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"


for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	input_dir=${root}/output_data/05_Taxanomic_Assignment
	output_dir=${root}/output_data/06_Generate_Output
	table_dir=${root}/output_data/03_Clustered_Data


	#Transform VSEARCH results to TSV
	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime tools export --input-path ${input_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.qza --output-path ${output_dir}
	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime tools export --input-path ${input_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.qza --output-path ${output_dir}
	echo mv ${output_dir}/taxonomy.tsv ${output_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv
	mv ${output_dir}/taxonomy.tsv ${output_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv

	#Transform nbClassifier results to TSV
	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime tools export --input-path ${input_dir}/${primer}_all_p985_taxa_nbclassified.qza --output-path ${output_dir}
	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime tools export --input-path ${input_dir}/${primer}_all_p985_taxa_nbclassified.qza --output-path ${output_dir}
	echo mv ${output_dir}/taxonomy.tsv ${output_dir}/${primer}_all_p985_taxa_nbclassified.tsv
	mv ${output_dir}/taxonomy.tsv ${output_dir}/${primer}_all_p985_taxa_nbclassified.tsv

	#Identify the taxa best classified by VSEARCH and NBClassifier
	echo singularity exec ${rimage} Rscript ${root}/code/06_Generate_Output/filterVsearch_nbClassifier.R ${output_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv ${output_dir}/${primer}_all_p985_taxa_nbclassified.tsv ${output_dir} ${primer}
	singularity exec ${rimage} Rscript ${root}/code/06_Generate_Output/filterVsearch_nbClassifier.R ${output_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.tsv ${output_dir}/${primer}_all_p985_taxa_nbclassified.tsv ${output_dir} ${primer}

   #Filter the VSEARCH file to retain just the Feature IDs required:
   	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime rescript filter-taxa \
	--i-taxonomy ${input_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.qza \
	--m-ids-to-keep-file ${output_dir}/${primer}_filtd_taxlist_vsearch.txt \
	--o-filtered-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_vsearch.qza
	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime rescript filter-taxa \
	--i-taxonomy ${input_dir}/${primer}_all_p985_taxa_VsearchOnly_p80_c94_COInr_Metazoa_and_Schmidt_LerayTrimmed.qza \
	--m-ids-to-keep-file ${output_dir}/${primer}_filtd_taxlist_vsearch.txt \
	--o-filtered-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_vsearch.qza

	#Filter the naive Bayes classifier samples to retain just the Feature IDs required:
	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime rescript filter-taxa \
	--i-taxonomy ${input_dir}/${primer}_all_p985_taxa_nbclassified.qza \
	--m-ids-to-keep-file ${output_dir}/${primer}_filtd_taxlist_sklearn.txt \
	--o-filtered-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_sklearn.qza
	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime rescript filter-taxa \
	--i-taxonomy ${input_dir}/${primer}_all_p985_taxa_nbclassified.qza \
	--m-ids-to-keep-file ${output_dir}/${primer}_filtd_taxlist_sklearn.txt \
	--o-filtered-taxonomy ${output_dir}/${primer}_all_p985_taxa_filtd_sklearn.qza

	#Combine the two taxonomy files and used the COMBINED lists of vsearch/sklearn-filtered data to filter the clustered seqs object
	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table merge-taxa \
	--i-data ${output_dir}/${primer}_all_p985_taxa_filtd_vsearch.qza ${output_dir}/${primer}_all_p985_taxa_filtd_sklearn.qza \
	--o-merged-data  ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza
	
	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table merge-taxa \
	--i-data ${output_dir}/${primer}_all_p985_taxa_filtd_vsearch.qza ${output_dir}/${primer}_all_p985_taxa_filtd_sklearn.qza \
	--o-merged-data  ${output_dir}/${primer}_all_p985_taxa_filtd_ALL.qza

	echo cut -f 1 -d ',' ${output_dir}/${primer}_filtd_tax_dataframe_ALL.csv > ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv
	cut -f 1 -d ',' ${output_dir}/${primer}_filtd_tax_dataframe_ALL.csv > ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv


	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table filter-seqs \
	--i-data ${table_dir}/${primer}_all_p985_seqs.qza \
	--o-filtered-data ${output_dir}/${primer}_all_p985_seqs_Filtd.qza \
	--m-metadata-file ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table filter-seqs \
	--i-data ${table_dir}/${primer}_all_p985_seqs.qza \
	--o-filtered-data ${output_dir}/${primer}_all_p985_seqs_Filtd.qza \
	--m-metadata-file ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv

	# Use the COMBINED lists of vsearch/sklearn-filtered data to filter the clustered table object
	echo singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table filter-features \
	--i-table ${table_dir}/${primer}_all_p985_table.qza \
	--o-filtered-table  ${output_dir}/${primer}_all_p985_table_filtd.qza \
	--m-metadata-file ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv

	singularity exec --bind $SCRATCH/tmp:/home/qiime2/q2cli $qiime_image qiime feature-table filter-features \
	--i-table ${table_dir}/${primer}_all_p985_table.qza \
	--o-filtered-table  ${output_dir}/${primer}_all_p985_table_filtd.qza \
	--m-metadata-file ${output_dir}/${primer}_filtd_tax_taxlist_ALL.csv
done