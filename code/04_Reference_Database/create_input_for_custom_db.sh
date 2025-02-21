#!/bin/bash
# Usage: ./script.sh <taxonomy_file> <fasta_file> [output_file]
# Example: ./script.sh tax_table.txt input.fasta output.tsv

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <taxonomy_file> <fasta_file> [output_file]"
  exit 1
fi

taxonomy_file="$1"
fasta_file="$2"
output_file="${3:-output.tsv}"

# Build an associative array from the taxonomy table.
# Key: first field (trimmed). Value: the rightmost nonblank field (trimmed).
declare -A tax_map

while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip empty lines.
  [[ -z "$line" ]] && continue

  # Trim leading/trailing whitespace.
  trimmed=$(echo "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')

  # Skip header line if first field contains "Vial" (or "Vial#")
  first_field=$(echo "$trimmed" | awk '{print $1}')
  if [[ "$first_field" == "Vial" || "$first_field" == "Vial#" ]]; then
    continue
  fi

  # Trim the key
  key=$(echo "$first_field" | xargs)

  # Use awk to iterate backwards over fields to find the first nonblank value,
  # and trim that field.
  tax=$(echo "$trimmed" | awk '{
      for(i = NF; i >= 1; i--) {
         gsub(/^[ \t]+|[ \t]+$/, "", $i);
         if($i != "") { print $i; exit }
      }
  }')
  tax_map["$key"]="$tax"
done < "$taxonomy_file"

# Process the FASTA file.
# For each FASTA record, we output:
#   Column 1: Full header (with comments removed; i.e., up to the first semicolon)
#   Column 2: Taxid (looked up via the key extracted from header; text before the first underscore)
#   Column 3: The concatenated sequence

# Initialize output file.
> "$output_file"

current_header=""
current_seq=""

while IFS= read -r line || [[ -n "$line" ]]; do
  # Remove any carriage return characters.
  line="${line//$'\r'/}"
  
  if [[ "$line" =~ ^\> ]]; then
    # If there is an existing record, output it.
    if [[ -n "$current_header" ]]; then
      # Remove the leading ">".
      full_header="${current_header#>}"
      # Remove the comment part (everything from the first semicolon onward).
      header_no_comment="${full_header%%;*}"
      # Extract the lookup key: text before the first underscore.
      key=$(echo "$header_no_comment" | cut -d'_' -f1 | xargs)
      taxid="${tax_map[$key]}"
      [[ -z "$taxid" ]] && taxid="NA"
      echo -e "${header_no_comment}\t${taxid}\t${current_seq}" >> "$output_file"
    fi
    # Start a new record.
    current_header="$line"
    current_seq=""
  else
    # Append sequence lines.
    current_seq+="$line"
  fi
done < "$fasta_file"

# Output the last record.
if [[ -n "$current_header" ]]; then
  full_header="${current_header#>}"
  header_no_comment="${full_header%%;*}"
  key=$(echo "$header_no_comment" | cut -d'_' -f1 | xargs)
  taxid="${tax_map[$key]}"
  [[ -z "$taxid" ]] && taxid="NA"
  echo -e "${header_no_comment}\t${taxid}\t${current_seq}" >> "$output_file"
fi

echo "Output written to $output_file"
