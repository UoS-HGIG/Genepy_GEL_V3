/*
    Module: Reatt_Genes
    Purpose:
        Reassemble per-gene `.meta` files that are split across aggV3 shard or
        subshard boundaries, so that each gene is scored once as a complete unit.

    Inputs:
        - CADD threshold identifier (`cadd`, currently 15 or 20)
        - Chromosome name (`chr`)
        - One or more folders containing per-gene `.meta` files

    Outputs:
        - Folder containing reassembled `.meta` files for genes found in more
          than one input folder (dup<cadd>)
        - List of genes that required reassembly (<chr>_<cadd>_dup.lst)

    Key steps:
        1. Create an output folder named `dup<cadd>` for the selected CADD
           threshold.
        2. Scan every input folder for `.meta` files.
        3. Group `.meta` files by gene name.
        4. For genes present in more than one folder, write a single combined
           `.meta` file:
           - Copy the header from the first file.
           - Append the variant rows from all files.
        5. Record each reassembled gene in the duplicate-gene list.

    Notes:
        - A gene may span a shard boundary, meaning different parts of the gene
          are processed in different subshards.
        - Reassembly prevents such genes from being scored more than once or
          from being scored using incomplete variant sets.
        - Genes found in only one input folder are not duplicated and are
          processed directly by the GenePy module.
*/
process Reatt_Genes {
    publishDir "${params.outDir}/${chr}", mode: "copy", overwrite: true
    //maxForks 10
    tag "Reatt_Genes"
    
    label "Reatt_Genes"
    input:
    tuple val(cadd),val(chr),path(folder_paths)
    //tuple path("metafilesALL"),path("metafiles15"),val("metafiles20")
    output:
    path(folder_paths), emit: path
    path("dup${cadd}"), emit: dup
    
    shell:
    """
    echo "start"
    
    OUTPUT_FOLDER="dup${cadd}"
    mkdir -p "\${OUTPUT_FOLDER}"
    touch "\${OUTPUT_FOLDER}/.keep"
    duplicated_genes="${chr}_${cadd}_dup.lst"
    > "\$duplicated_genes"
    set -u
    declare -a FOLDERS
 

    declare -A gene_files
    for dir in ${folder_paths}; do 
        echo "\$dir"
        clean_paths=\$(echo \$dir | tr -d '[],')
        echo "\$clean_paths"
        ls \$clean_paths
        if compgen -G "\$clean_paths"/*.meta > /dev/null; then
            for file in "\$clean_paths"/*.meta; do
              ##  echo "File: \$file" >> "\$GENE_LIST"
                gene_name=\$(basename "\$file")
                gene_name="\${gene_name%.meta}"
                gene_files["\$gene_name"]="\${gene_files[\$gene_name]:-} \$file"
                echo "\$gene_name"
            done
        else
            echo "No .meta files found in \$clean_paths"
        fi
    done
    

    echo "half done!"
    for item in "\${gene_files[@]}"; do
        echo "\$item"
    done

    for gene in "\${!gene_files[@]}"; do
        
        echo "this is  \$gene"
        files=(\${gene_files[\$gene]})
        echo "\$files"
        if [ \${#files[@]} -gt 1 ]; then
            echo "this is duplicated  \$gene"
            echo "\$gene" >> "\$duplicated_genes"
            output_file="\${OUTPUT_FOLDER}/\${gene}.meta"

            head -n 1 "\${files[0]}" > "\$output_file"

            for file in "\${files[@]}"; do
                tail -n +2 "\$file" >> "\$output_file"
            done
            ##echo "\$output_file" >> "\$OUTPUT_FILE_LIST" 
        fi
    done
   
    """
}
