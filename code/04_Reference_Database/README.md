# Download COInr Database

	source gut_contents.env
	cd ${root}/input_data/04_Create_Reference_Database
	wget https://zenodo.org/records/6555985/files/COInr_2022_05_06.tar.gz?download=1
	tar -zxvf COInr_2022_04_06.tar.gz
	mv ${root}/input_data/04_Reference_Database/COInr_2022_05_06/COInr.tsv ${root}/input_data/04_Create_Reference_Database/
	mv ${root}/input_data/04_Reference_Database/taxonomy.tsv ${root}/input_data/04_Create_Reference_Database/COInr.taxonomy.tsv
	rm -rf ${root}/input_data/04_Reference_Database/COInr_2022_05_06

# Create A Filtered COI+NCBI NR Database for Metazoans with a Minimum Taxanomic Level of Family

	source gut_contents.env

	export input_data=${root}/input_data/04_Reference_Database/
	export output_data=${root}/output_data/04_Reference_Database/

	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/select_taxa.pl -taxon_list ${input_data}/taxa_whitelist4mkCOInr.txt -tsv ${input_data}/COInr.tsv -taxonomy ${input_data}/COInr.taxonomy.tsv -outdir ${output_data} -out COInr.taxfilt.tsv -min_taxlevel family


# Process Schmidt's Sequenced Barcodes

Total 358 barcodes passed QC. For Run2, Ray only sent me all the barcodes and not just the ones that passed (175 total, 104 passed). So I am adding 429 barcodes total.

	source gut_contents.env

	export input_data=${root}/input_data/04_Reference_Database/schmidt_barcodes/Sequence_Files/
	export output_data=${root}/output_data/04_Reference_Database/
	export taxref=${root}/input_data/04_Reference_Database/schmidt_barcodes/Gabon_inverts_taxonomy_ERIC_JG_CORRECTED.txt

	#Change mutliline run2 fasta file to single line.
	awk '{gsub(/\r/,"");} 
     /^>/{if(seq) print seq; print; seq=""; next} 
     {seq = seq $0} 
     END{print seq}' ${input_data}/GabonRun2_consensus_all_step1.fas > ${input_data}/GabonRun2_consensus_all_step1.fa


	# Append "_2" to the name of any duplicates

		awk '{
		gsub(/\r/,"")
	}
	/^>/{
		split($0, a, ";")
		key = a[1]
		if (seen[key]++){
			newheader = a[1] "_2"
			for(i=2; i<=length(a); i++){
				newheader = newheader ";" a[i]
			}
			print newheader
		} else {
			print $0
		}
		next
	}
	{ print }' "${input_data}/GabonRun2_consensus_all_step1.fa" > "${input_data}/GabonRun2_consensus_all_step1.duprename.fa"


	#Merge Everything
	cat ${input_data}/*predgood_barcodes.fa ${input_data}/GabonRun2_consensus_all_step1.duprename.fa > ${input_data}/allschmidt429.fasta

	#Check that we have 429 sequences
 	grep ">" ${input_data}/allschmidt429.fasta | wc -l 

	#NB: I think Ray and the lady may have gotten it wrong and wrote 3212 instead of 3612. I am physically changing 3212 to 3612 in the allschmidt429.fasta file. 
	sed -i 's/3212/3612/g' ${input_data}/allschmidt429.fasta


	#Convert to TSV
	bash ${root}/code/04_Reference_Database/create_input_for_custom_db.sh ${taxref} ${input_data}/allschmidt429.fasta ${output_data}/allschmidt429.tsv

	#Process Schmidt Barcodes with COInr
	export input_data=${root}/input_data/04_Reference_Database/
	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/format_custom.pl -custom ${output_data}/allschmidt429.tsv -taxonomy ${input_data}/COInr.taxonomy.tsv -outdir ${output_data}

	#Manually remove 'homonymy'=1 rows that have lower phylogenetic resolution than others, and delete the column for downstream processing.  save as "custom_lineages_verified.tsv"

	#Add TaxIDs to Custom Sequences

	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/add_taxids.pl -lineages ${output_data}/custom_lineages_verified.tsv -sequences ${output_data}/custom_sequences.tsv -taxonomy ${input_data}/COInr.taxonomy.tsv -outdir ${output_data}

	# Dereplicate Custom Sequences

	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/dereplicate.pl -tsv ${output_data}/sequences_with_taxIDs.tsv -outdir ${output_data} -out custom_dereplicated_sequences.tsv


## Merge Schmidt and COInr Metazoa Database

	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/pool_and_dereplicate.pl -tsv1 ${output_data}/COInr.taxfilt.tsv -tsv2 ${output_data}/custom_dereplicated_sequences.tsv -outdir ${output_data} -out COInr_custom.tsv


## Restrict to the Leray Region
	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/select_region.pl -tsv ${output_data}/COInr_custom.tsv -outdir ${output_data} -e_pcr 1 -fw GCHCCHGAYATRGCHTTYCC -rv TCDGGRTGNCCRAARAAYCA -trim_error 0.3 -min_amplicon_length 280 -max_amplicon_length 345 -min_overlap 20 -tcov 0.9 -identity 0.7

## Make Qiime Artifact
	singularity exec ${mkcoinr_image} perl ${mkcoinr_path}/format_db.pl -tsv ${output_data}/trimmed.tsv -taxonomy ${output_data}/taxonomy_updated.tsv -outfmt qiime -outdir ${output_data} -out COInr_Metazoa_and_Schmidt_LerayTrimmed.qiime

## Run Training
	bash code/04_Reference_Database/01_run_classifier_training.sh