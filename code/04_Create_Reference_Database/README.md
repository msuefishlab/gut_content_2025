
JG: This is sort of a mess.  Going to need the text file from sophie called:

"Gabon_inverts_taxonomy_ERIC.txt"

Going to trust this one for the time being, but will ultimately need to reexamine this.




3. Creating the reference database: Merging O’rouke et al. 2021 bigCOI and the newly barcoded bugs from Schmidt 
1. Download the bigCOI database from O’Rourke

“We collected reference sequences and associated taxonomy information from two resources: BOLD (Ratnasingham and Hebert, 2007) and a GenBank-derived dataset curated by Terri Porter (Porter and Hajibabaei, 2018). Reference sequences included COI records from arthropod, chordate, and other animal taxa, as well as fungal, protist, and other microeukaryote COI records. We dereplicated the initial collection of sequences, then applied a Least Common Ancestor
(LCA) process using a consensus approach to classify records that shared identical sequence information but differed with respect to taxonomic information. Additional filters included discarding references with non-standard IUPAC DNA characters, removing sequences <100 bp, and retaining only references that contained at least family-level names. The final dataset included 2,181,331 distinct sequences. The construction of this database is described here: https://github.com/devonorourke/mysosoup/blob/master/docs/database_construction.md.

Downloaded bigCOI.derep.seqs.qza and bigCOI.derep.tax.qza from 
https://osf.io/qju3w/

There is a lot of ambiguous taxa that has no resolution still. 480,655 sequences have this label:
Ambiguous_taxa;Ambiguous_taxa;c__Ambiguous_taxa;o__Ambiguous_taxa;f__Ambiguous_taxa;g__Ambiguous_taxa;s__Ambiguous_taxa 
So let’s remove it

    ## filter taxonomy file
    qiime rescript filter-taxa \
    --i-taxonomy bigCOI.derep.tax.qza \
    --p-exclude 'Ambiguous_taxa;Ambiguous_taxa;c__Ambiguous_taxa;o__Ambiguous_taxa;f__Ambiguous_taxa;g__Ambiguous_taxa;s__Ambiguous_taxa' \
    --o-filtered-taxonomy bigCOI.derep.tax.Filtd1.qza
    
    ## filter sequences files so that the ones that are Ambiguous are out (n=480,655!)
    qiime taxa filter-seqs \
    --i-sequences bigCOI.derep.seqs.qza \
    --i-taxonomy bigCOI.derep.tax.qza \
    --p-exclude 'Ambiguous_taxa;Ambiguous_taxa;c__Ambiguous_taxa;o__Ambiguous_taxa;f__Ambiguous_taxa;g__Ambiguous_taxa;s__Ambiguous_taxa' \
    --o-filtered-sequences bigCOI.derep.seqs.Filtd1.qza
    
    qiime metadata tabulate  \
    --m-input-file  bigCOI.derep.tax.Filtd1.qza \
    --o-visualization bigCOI.derep.tax.Filtd1.qzv
    
    #Get at least class information (n=458)
    qiime rescript filter-taxa \
    --i-taxonomy bigCOI.derep.tax.Filtd1.qza \
    --p-exclude 'c__c__;o__o__;f__f__;g__g__;s__s__' \
    --o-filtered-taxonomy bigCOI.derep.tax.Filtd2.qza
    
    ## filter sequences file
    qiime taxa filter-seqs \
    --i-sequences  bigCOI.derep.seqs.Filtd1.qza \
    --i-taxonomy bigCOI.derep.tax.Filtd1.qza \
    --p-exclude 'c__c__;o__o__;f__f__;g__g__;s__s__' \
    --o-filtered-sequences bigCOI.derep.seqs.Filtd2.qza
    
    qiime metadata tabulate  \
    --m-input-file  bigCOI.derep.tax.Filtd2.qza \
    --o-visualization bigCOI.derep.tax.Filtd2.qzv
    # The Filtd2 database has 1,600,953 sequences. It’s super huge and takes >48 hours and 256GB Ram to train  Bayesian classifier.
    
    # Instead, we are going to go with this one that has family level information. We are also going to trim to primer regions. 
    
    #Get at least family level information (n=67135)
    qiime rescript filter-taxa \
    --i-taxonomy bigCOI.derep.tax.Filtd1.qza \
    --p-exclude 'f__f__;g__g__;s__s__' \
    --o-filtered-taxonomy bigCOI.derep.tax.FiltdFamily.qza
    
    # filter sequences file
    qiime taxa filter-seqs \
    --i-sequences  bigCOI.derep.seqs.Filtd1.qza \
    --i-taxonomy bigCOI.derep.tax.Filtd1.qza \
    --p-exclude 'f__f__;g__g__;s__s__' \
    --o-filtered-sequences bigCOI.derep.seqs.FiltdFamily.qza
    
    qiime metadata tabulate  \
    --m-input-file  bigCOI.derep.tax.FiltdFamily.qza \
    --o-visualization bigCOI.derep.tax.FiltdFamily.qzv
    #1,534,276 total sequences

 


2. Get all the data back from Ray Schmidt

Total 358 barcodes passed QC. For Run2, Ray only sent me all the barcodes and not just the ones that passed (175 total, 104 passed). So I am adding 429 barcodes total.

1.Merge all of them into one document.

First change interleaved run2  fasta file to single line. This was really hard, I ended up downloading seqinr package in r and running this

    install.packages('seqinr') 
    library(seqinr)
    seqs = read.fasta(file='GabonRun2_consensus_all_step1.fas')
    write.fasta(seqs, names(seqs), nbchar=658, file.out='sequences2.fa')
    
    ## then run this in unix
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' sequences2.fa > sequences3.fa 
    
    #renamed sequences3 to GabonRun2_consensus_all_step1_oneline.fas

There are 21 repeat sequence names in the GabonRun2 file. I add a _2 to the repeats. No idea why thats there (did it in the oneline file)

    3592_13_all
    3592_14_all
    3592_15_all
    3592_16_all
    3592_17_all
    3592_18_all
    3592_23_all
    3592_31_all
    3592_35_all
    3592_37_all
    3592_38_all
    3592_40_all
    3592_41_all
    3592_43_all removed the second one it was full of NNN
    3592_46_all
    3592_4_all
    3592_6_all
    3592_7_all
    3592_8_all
    3598_1_all
    3598_2_all

 Navigate to the folder where we have the four different fasta files from the four different runs and run this:

    cat * > allschmidt429.fasta
    #Check that we have 429 sequences
     grep ">" allschmidt429.fasta | wc -l 

NB: I think Ray and the lady may have gotten it wrong and wrote 3212 instead of 3612. I am physically changing 3212 to 3612 in the allschmidt429.fasta file. 
2.Create taxonomy file and upload taxonomy and sequence files into QIIME2


    # Access the names of his sequences
    grep -o -E "^>\w+" allschmidt429.fasta | tr -d ">" >names.txt
    #Access just the numbers
    cut -d "_" -f1 names.txt >numbers.txt
    #access just the sequences
    cat allschmidt429.fasta |grep -v "^>" > seqsschmidt.tmp
    
    paste names.txt seqsschmidt.tmp | sed 's/^/>/' | tr '\t' '\n' > schmidt_all_seqs.fa
    # Paste the numbers as the column next to it (got rid of the quality stuff)
    paste -d"\t" names.txt numbers.txt > seqs_names.txt
    
    #get the taxonomy info for each sequence
    awk ' FILENAME=="Gabon_inverts_taxonomy_ERIC.txt" {arr[$1]=$0; next}
            FILENAME=="seqs_names.txt"  {print arr[$2]} ' Gabon_inverts_taxonomy_ERIC.txt  seqs_names.txt > schmidt_all_taxonomy.txt
    
    #paste all taxonomy info together
    awk 'BEGIN { IFS=OFS="\t" }
         NR >0 { two=$2
                 for(i=3; i <= NF; i++) two=two";"$i
                 print $1 OFS two }' schmidt_all_taxonomy.txt > schmidt_taxmap_all.txt
    
    #add each seqs name
    #first, sort
    sort -k1 schmidt_taxmap_all.txt > sorted_schmidt_taxmap_all.txt
    sort -k 2 seqs_names.txt   > sorted_seqs_names.txt
    #then did this manually in excel because UGH final taxamp file called schmidt_taxmap_all_final.txt
    #With Cyberduck, add the fasta and taxonomy files schmidt_all_seqs.fa and schmidt_may2_taxmap_final.txt 

Doing this in R instead in Schmidt_taxonomy.Rmd

3.First convert the schmidt sequences to qza

    conda activate qiime2-2022.2
    ## fasta file import
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path schmidt_all_429seqs.fa \
      --output-path schmidt_all_429seqs.qza


3. Merge with the bigCOI sequences and same for taxonomies
    qiime feature-table merge-seqs \
      --i-data bigCOI.derep.seqs.FiltdFamily.qza \
      --i-data ~/Schmidt_gabon_COI/schmidt_all_429seqs.qza \
      --o-merged-data bigCOI_FiltdFamily_Schmidtfinal428.qza
    
    ## taxonomy file import
     qiime tools import \
       --type 'FeatureData[Taxonomy]' \
       --input-format HeaderlessTSVTaxonomyFormat \
      --input-path Schmidt428_taxo_for_qiime.txt \
       --output-path Schmidt428_taxo_for_qiime.qza
    
    ## Merge with with the bigCOI taxonomy
     qiime feature-table merge-taxa \
      --i-data ~/Schmidt_gabon_COI/Schmidt428_taxo_for_qiime.qza \
      --i-data bigCOI.derep.tax.FiltdFamily.qza \
      --o-merged-data bigCOI_SChmidt428_taxonomy_Familyfinal.qza
    
    qiime metadata tabulate  \
    --m-input-file  bigCOI_FiltdFamily_Schmidtfinal428.qza \
    --m-input-file  bigCOI_SChmidt428_taxonomy_Familyfinal.qza \
    --o-visualization  bigCOIFamily_SChmidt428_seqtaxonomy_final.qzv
    


4. Extract primer regions
    #Primerset2
    qiime feature-classifier extract-reads \
       --i-sequences bigCOI_FiltdFamily_Schmidtfinal428.qza \
       --p-f-primer GCHCCHGAYATRGCHTTYCC \
       --p-r-primer ARYATDGTRATDGCHCCDGC \
       --p-n-jobs 8 \
       --o-reads bigCOI_FiltdFamily_Schmidtfinal428_primerset2_seq.qza
    
    #Primerset1
    qiime feature-classifier extract-reads \
       --i-sequences bigCOI_FiltdFamily_Schmidtfinal428.qza \
       --p-f-primer GGWACWGGWTGAACWGTWTAYCCYCC \
       --p-r-primer TCDGGRTGNCCRAARAAYCA \
       --p-n-jobs 8 \
       --o-reads bigCOI_FiltdFamily_Schmidtfinal428_primerset1_seq.qza


5. Dereplicate again
    #Primerset2
    qiime rescript dereplicate \
        --i-sequences bigCOI_FiltdFamily_Schmidtfinal428_primerset2_seq.qza \
        --i-taxa bigCOI_SChmidt428_taxonomy_Familyfinal.qza \
        --p-mode 'lca' \
        --o-dereplicated-sequences bigCOI_FiltdFamily_Schmidtfinal428_primerset2_derep.qza \
        --o-dereplicated-taxa  bigCOI_SChmidt428_taxonomy_Familyfinal_primerset2_derep.qza
    
    qiime tools export --input-path bigCOI_FiltdFamily_Schmidtfinal428_primerset2_derep.qza --output-path exportation/
    
    mv dna-sequences.fasta bigCOI_FiltdFamily_Schmidtfinal428_primerset2_derep.fasta
    #How many sequences
    grep -c "^>" bigCOI_FiltdFamily_Schmidtfinal428_primerset1_derep.fasta
    #821327
    
    #Primerset1
    qiime rescript dereplicate \
        --i-sequences bigCOI_FiltdFamily_Schmidtfinal428_primerset1_seq.qza \
        --i-taxa bigCOI_SChmidt428_taxonomy_Familyfinal.qza \
        --p-mode 'lca' \
        --o-dereplicated-sequences bigCOI_FiltdFamily_Schmidtfinal428_primerset1_derep.qza \
        --o-dereplicated-taxa  bigCOI_SChmidt428_taxonomy_Familyfinal_primerset1_derep.qza
    
    qiime tools export --input-path bigCOI_FiltdFamily_Schmidtfinal428_primerset1_derep.qza --output-path exportation/
    
    mv dna-sequences.fasta bigCOI_FiltdFamily_Schmidtfinal428_primerset1_derep.fasta
    grep -c "^>" bigCOI_FiltdFamily_Schmidtfinal428_primerset1_derep.fasta
    #364321
