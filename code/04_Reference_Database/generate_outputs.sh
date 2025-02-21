#!/bin/bash
# Usage: ./generate_outputs.sh <filtered_tsv_file> <output_base_name>
# This script produces four output files based on a user-specified base name:
#   <base>.seqs.csv        - CSV file with processid (as sequenceID) and nucleotides
#   <base>.taxa_string.csv   - CSV file with processid and a taxonomy string
#   <base>.meta.txt          - Semicolon-delimited metadata file with quoted fields
#   <base>.fasta             - FASTA formatted file with headers ">seqID tax=taxon" followed by nucleotides

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <filtered_tsv_file> <output_base_name>"
    exit 1
fi

input="$1"
base="$2"

awk -F'\t' -v base="$base" 'BEGIN {
    # Define output file names based on the provided base name.
    seq_out   = base ".seqs.csv";
    tax_out   = base ".taxa_string.csv";
    meta_out  = base ".meta.txt";
    fasta_out = base ".fasta";

    # Write headers to the CSV and metadata output files.
    print "sequenceID,nucleotides" > seq_out;
    print "sequenceID,taxon" > tax_out;
    print "\"sequenceID\";\"processid\";\"bin_uri\";\"genbank_accession\";\"country\";\"institution_storing\";\"phylum\";\"class\";\"order\";\"family\";\"genus\";\"species\"" > meta_out;
}

NR==1 {
    # Skip the header line of the input TSV.
    next;
}

{
    # Map columns based on the BOLD header:
    # processid (sequenceID)          = $1   (used as sequenceID)
    # nuc (nucleotides)               = $67
    # kingdom                         = $14
    # phylum                          = $15
    # class                           = $16
    # order                           = $17
    # family                          = $18
    # genus                           = $21
    # species                         = $22
    # bin_uri                         = $8
    # insdc_acs (genbank_accession)   = $69
    # country/ocean                   = $48
    # inst (institution_storing)      = $11

    seqID       = $1;
    nucleotides = $67;

    # Build the taxon string with prefixes.
    kingdom = "k__" $14;
    phylum  = "p__" $15;
    class   = "c__" $16;
    order   = "o__" $17;
    family  = "f__" $18;
    genus   = "g__" $21;
    species = "s__" $22;
    taxon   = kingdom ";" phylum ";" class ";" order ";" family ";" genus ";" species;

    # Write to the sequences CSV: sequenceID,nucleotides.
    print seqID "," nucleotides > seq_out;

    # Write to the taxonomy CSV: sequenceID,taxon.
    print seqID "," taxon > tax_out;

    # Write to the metadata file (semicolon-delimited with each field quoted).
    # Using processid ($1) as both sequenceID and processid.
    printf("\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\";\"%s\"\n", \
        seqID, $1, $8, $69, $48, $11, $15, $16, $17, $18, $21, $22) > meta_out;

    # Write to the FASTA file: header as ">seqID tax=taxon" followed by nucleotides.
    print ">" seqID ";tax=" taxon > fasta_out;
    print nucleotides > fasta_out;
}
END {
    close(seq_out);
    close(tax_out);
    close(meta_out);
    close(fasta_out);
}' "$input"
