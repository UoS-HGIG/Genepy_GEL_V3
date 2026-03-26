input_file= $1
header_file="header.txt"
head -n 1 "$input_file" > "$header_file"
mkdir gene_files

# split by samples
tot_c=$(awk '{print NF - 26; exit}' $header_file)
gs=$((tot_c / 3))

mkdir gp1
mkdir gp2
mkdir gp3

cut -d$'\t' -f1-26,27-$((26+gs)) $header_file > gp1/header.txt
cut -d$'\t' -f1-26,$((26+gs))-$((26+ 2*gs)) $header_file > gp2/header.txt
cut -d$'\t' -f1-26,$((26+ 2*gs))- $header_file > gp3/header.txt


awk -F'\t' -v header="$(cat "$header_file")" 'NR > 1 {
    gene_name = $6
    gene_info = $0  # Entire line is gene info
    if (match(gene_name, /ENSG[0-9]+(\.[0-9]+)?/)) {
        gene_id = substr(gene_name, RSTART, RLENGTH)
        gene_file = "gene_files/" gene_id ".txt"
        if (!seen[gene_file]++) {
            # Print header to gene file
            print header > gene_file
        }
        print gene_info >> gene_file
    }
}' "$input_file"

gene_files_folder="gene_files"
 ls "$gene_files_folder"/*.txt | sed 's|.*/\(.*\)\.txt|\1|' > gene_list.txt
 input_gene_list="gene_list.txt"
 output_gene_size="gene_sizes.txt"
 > "$output_gene_size" # Clear the output file
while IFS= read -r gene; do     cadd_file="gene_files/${gene}.txt";     if [ ! -f "$cadd_file" ]; then         echo "Warning: ${gene}.txt file not found for gene '$gene'.";         continue;     fi;      size=$(stat -c "%s" "$cadd_file");      echo "$size $gene" >> "$output_gene_size"; done < "$input_gene_list"
sort -n -o "$output_gene_size" "$output_gene_size"

cat gene_list.txt | while read i; do
    cp gp2/header.txt gp2/${i}.meta; cp gp1/header.txt gp1/${i}.meta; cp gp3/header.txt gp3/${i}.meta; cat gene_files/${i}.txt | cut -d$'\t' -f1-26,27-$((26+gs)) > gp1/${i}.meta; cat gene_files/${i}.txt | cut -d$'\t' -f1-26,$((26+gs))-$((26+ 2*gs)) > gp2/${i}.meta; cat gene_files/${i}.txt | cut -d$'\t' -f1-26,$((26+ 2*gs))- > gp3/${i}.meta;
done