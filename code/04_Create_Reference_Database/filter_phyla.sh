#!/bin/bash
# Usage: ./filter_phyla.sh whitelist_file.tsv BOLD_database.tsv filtered_output.tsv

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <whitelist_file> <tsv_file> <output_file>"
    exit 1
fi

whitelist_file="$1"
tsv_file="$2"
output_file="$3"

awk -F'\t' '
NR==FNR {
    # Build the whitelist from the first file (one taxon per line)
    whitelist[$1]=1;
    next
}
FNR==1 {
    # Always print the header from the TSV file
    print;
    next
}
# For other lines, print only if:
#   - The phylum (column 15) is in the whitelist AND
#   - The marker_code (column 71) equals "COI-5P"
(($15 in whitelist) && ($71 == "COI-5P"))
' "$whitelist_file" "$tsv_file" > "$output_file"
