#!/bin/bash
# Usage: ./filter_sequences.sh input.tsv output.tsv
# The input TSV is assumed to be the output from the previous filtering step.
# This script applies additional sequence and taxonomy filters.

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_tsv> <output_tsv>"
    exit 1
fi

input="$1"
output="$2"

awk -F'\t' 'BEGIN { OFS="\t" }
NR==1 {
    # Print header unchanged
    print;
    next;
}
{
    # Make a copy of the nuc sequence (column 67)
    seq = $67;
    # Filter 1: Remove gap characters ("-") and check sequence length (> 100bp)
    gsub(/-/, "", seq);
    if (length(seq) <= 100)
        next;
    
    # Filter 2: Remove sequences with any characters other than A, T, G, or C
    if (seq !~ /^[ATGC]+$/)
        next;
    
    # Filter 3: For Arthropoda or Chordata, require the family field (column 18) be non-empty.
    if (($15 == "Arthropoda" || $15 == "Chordata") && ($18 == ""))
        next;
    
    # Update the nuc field with the gap-free sequence
    $67 = seq;
    
    # Print the modified record
    print;
}' "$input" > "$output"
