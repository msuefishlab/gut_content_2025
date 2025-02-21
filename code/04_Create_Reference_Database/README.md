#Porter Lab v 5 COI Database:

    wget https://github.com/terrimporter/CO1Classifier/releases/download/RDP-COI-v5.0.0-ref/RDP_COIv5.0.0_ref.zip

    unzip RDP_COIv5.0.0_ref.zip

    cd mydata_ref

    cat mytrainseq.fasta | paste - - | grep -v "Arthropoda" | grep -v "Chordata" | tr '\t' '\n' > tPorter.nonArth_nonChord.fasta

    cat tPorter.nonArth_nonChord.fasta | paste - - | cut -f 1 -d ' ' | sed 's/>//' > left.tmp
    cat tPorter.nonArth_nonChord.fasta | grep '^>' | cut -f 2- -d ' ' | cut -f 3- -d ';' > right.tmp
    paste left.tmp right.tmp > tPorter.nonArth_nonChord.taxa

    cat tPorter.nonArth_nonChord.fasta | grep -v "^>" > bottom.tmp
    paste left.tmp bottom.tmp | sed 's/^/>/' | tr '\t' '\n' > tPorter.nonArth_nonChord.nolabel.fasta

    mv tPorter.nonArth_nonChord.taxa ${root}/output_data/04_Reference_Database
    mv tPorter.nonArth_nonChord.nolabel.fasta ${root}/output_data/04_Reference_Database


# BOLD

Download BOLD Database (24-Jan-2025) and upload to ${root/input_data/04_Create_Reference_Database}

    source gut_contents.env 

## Filter for Taxa:
    bash ${root}/code/04_Create_Reference_Database/filter_phyla.sh ${root}/input_data/04_Reference_Database/taxa_whitelist.txt ${root}/input_data/04_Reference_Database/BOLD_Public.24-Jan-2025.fasta  ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.tsv

## Filter Sequences:
    bash ${root}/code/04_Create_Reference_Database/filter_sequences.sh ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.tsv ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.tsv

## Create Output Files:
    bash ${root}/code/04_Create_Reference_Database/generate_outputs.sh ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.tsv ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.formatted

## Taxonomy mapping file:
    cat ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.formatted.fasta | grep '^>' | sed 's/^>//' | cut -d ';' -f 1 > ${root}/output_data/04_Reference_Database/tmp.left
    cat ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.formatted.fasta  | grep '^>' | sed 's/^>//' | cut -d ';' -f2- | sed 's/tax=k__//' | sed 's/p__//' | sed 's/c__//' | sed 's/o__//' | sed 's/f__//' | sed 's/g__//' | sed 's/s__//' > ${root}/output_data/04_Reference_Database/tmp.right
    paste -d '\t'  ${root}/output_data/04_Reference_Database/tmp.left ${root}/output_data/04_Reference_Database/tmp.right | gzip > ${root}/output_data/04_Reference_Database/tmp_nolabels.taxa.gz
    rm ${root}/output_data/04_Reference_Database/tmp.left ${root}/output_data/04_Reference_Database/tmp.right

## reduced fasta file:
    cat ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.formatted.fasta | grep '^>' | cut -d ';' -f 1 > ${root}/output_data/04_Reference_Database/tmp.left
    cat ${root}/output_data/04_Reference_Database/BOLD_Public.24-Jan-2025.taxfilt.seqfilt.formatted.fasta | grep -v '^>' > ${root}/output_data/04_Reference_Database/tmp.seqs
    paste -d '\t' ${root}/output_data/04_Reference_Database/tmp.left ${root}/output_data/04_Reference_Database/tmp.seqs | tr '\t' '\n' | gzip > ${root}/output_data/04_Reference_Database/tmp_nolabels.fasta.gz
    rm ${root}/output_data/04_Reference_Database/tmp.left ${root}/output_data/04_Reference_Database/tmp.seqs



# Combine Everything:

    cat ${root}/output_data/04_Reference_Database/tPorter.nonArth_nonChord.taxa.gz ${root}/output_data/04_Reference_Database/tmp_nolabels.taxa.gz > allRecords_taxmap.txt.gz
    cat ${root}/output_data/04_Reference_Database/tPorter.nonArth_nonChord.nolabel.fasta.gz ${root}/output_data/04_Reference_Database/tmp_nolabels.fasta.gz > allRecords_nolabels.fasta.gz


# Deduplicate

