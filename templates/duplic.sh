#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 folder_list.txt chromosome_name output_path" 
    exit 1
fi

folder_list="$1"
chromosome_name="$2"
output_path="$3"
# Create a directory for concatenated genes
OUTPUT_FOLDER="$output_path/concatenated_genes_${chromosome_name}/metafiles"
mkdir -p "$OUTPUT_FOLDER"

FINAL_LIST="${output_path}/${chromosome_name}_final.lst"


# Read folder list into an array
mapfile -t FOLDERS < "$folder_list"

# Declare an associative array to track duplicate genes
declare -A gene_files
> "$FINAL_LIST"
# Loop through each folder and find duplicate genes
for folder in "${FOLDERS[@]}"; do
echo "$folder"/metafiles/metafiles >> "$FINAL_LIST"
    for file in "$folder"/metafiles/metafiles/*.meta; do
        [ -e "$file" ] || continue  # Skip if no .meta files found
        gene_name=$(basename "$file")
        gene_files["$gene_name"]+="$file "
    done
    
done

# Loop through the associative array and concatenate files for duplicate genes
for gene in "${!gene_files[@]}"; do
    files=(${gene_files[$gene]})
    if [ ${#files[@]} -gt 1 ]; then
        output_file="$OUTPUT_FOLDER/${gene}"
        
        # Extract the header from the first file
        head -n 1 "${files[0]}" > "$output_file"

        # Concatenate the rest of the files excluding their headers
        for file in "${files[@]}"; do
            tail -n +2 "$file" >> "$output_file"
        done
        
        ## Remove original files
        for file in "${files[@]}"; do
            rm -r "$file"
        done
    fi
done

if [ -d "$OUTPUT_FOLDER" ]; then
  echo "${output_path}/concatenated_genes_${chromosome_name}/metafiles" >> "$FINAL_LIST"
fi
