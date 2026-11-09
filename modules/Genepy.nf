/*
    Module: Genepy_score
    Purpose:
        Run GenePy to calculate a per-gene pathogenicity score from the
        per-gene `.meta` matrices produced by the pipeline.

    Inputs:
        - Folder(s) containing per-gene `.meta` matrices
        - Chromosome name (`chr`)
        - CADD threshold identifier (`cadd`); only the CADD15 folders are passed in,
          and genepy.py writes both CADD15 and CADD20 scores from them
        - GenePy Python script (`genepy.py`)
        - Karyotype/sex file (`kary`)
        - Folder of reassembled cross-shard genes (`dup`)

    Outputs:
        - Six GenePy output files per gene, each sample / score / label:
            `*_CADD15.meta.txt`, `*_CADD20.meta.txt`                (whole gene, as before)
            `*_CADD15_core.meta.txt`, `*_CADD20_core.meta.txt`      (variants inside the gene boundary)
            `*_CADD15_CRE.meta.txt`, `*_CADD20_CRE.meta.txt`        (variants within cre_flank of the most
                                                                     upstream TSS, outside the gene)
          published to Genepy_score/<cadd>/, Genepy_score/<cadd>/core/ and Genepy_score/<cadd>/CRE/

    Key steps:
        1. Copy the GenePy Python script into the task working directory.
        2. For every `.meta` file:
           - Replace empty or zero-like allele-frequency values in columns 7-16
             with a small minimum frequency value.
           - Run `genepy.py` with the processed `.meta` file and karyotype file.
        3. If a gene has a reassembled version in the `dup` folder, use that
           complete cross-shard file instead of the individual subshard file.

    Notes:
        - Columns 7-16 correspond to allele-frequency slots for up to ten ALT
          alleles.
        - The frequency-flooring step prevents zero or missing frequencies from
          causing problems in the GenePy calculation.
        - The `.meta` column layout must remain consistent with the output of
          `templates/pre_1.sh` / `pre_1g_v2.sh`; changing the number of allele
          frequency columns requires updating this module as well.
        - The final GenePy scores are published under
          `<outDir>/<chr>/Genepy_score/<cadd>/`.
*/
process Genepy_score {
    // One task scores every gene at CADD15 and CADD20 from the CADD15 .meta matrices
    // (CADD20 = same rows with raw CADD < 2.097252 blanked; done inside genepy.py).
    publishDir "${params.outDir}/${chr}/Genepy_score", mode: "copy", overwrite: true,
        saveAs: { f ->
            def m = (f =~ /_CADD(\d+)(?:_(core|CRE))?\.meta\.txt$/)
            m ? (m[0][2] ? "${m[0][1]}/${m[0][2]}/${f}" : "${m[0][1]}/${f}") : f
        }
    tag "Genepy_score_${chr}"
    label "Genepy_score"
    maxForks 20
    input:
    tuple path(path1),val(chr),val(cadd),path(genepy),path(kary),path(dup)
    path(gene_bed)   // BED6 gene boundaries with strand; core and CRE are derived from it

    output:
    path("*_CADD15.meta.txt"), optional: true, emit: cadd15
    path("*_CADD20.meta.txt"), optional: true, emit: cadd20
    path("*_CADD15_core.meta.txt"), optional: true, emit: cadd15_core
    path("*_CADD15_CRE.meta.txt"),  optional: true, emit: cadd15_cre
    path("*_CADD20_core.meta.txt"), optional: true, emit: cadd20_core
    path("*_CADD20_CRE.meta.txt"),  optional: true, emit: cadd20_cre

    shell:
    """
    Genepy=\$(readlink -f ${genepy})
    ls \$Genepy
    ls ${path1}
    cp \$Genepy ./gp.py
    chmod +x ./gp.py
   
    for file in ${path1}/*; do
            if [ -f "\$file" ]; then
                fname=\$(basename "\$file")
                if [[ ! "\$fname" == *.meta ]]; then
                    echo "Skipping \$fname (not a .meta file)"
                    continue
                fi
                # Case 1: path1 itself is a dup folder → always process
                if [ -n "${dup}" ] && [ -f "${dup}/\$fname" ]; then
                    echo "Processing \$file from dup folder"
                     awk -F"\\t" '{OFS=FS;for (i=7;i<=16;i++) { if(length(\$i)<1 || \$i ~ /^0+([.0]+)?([eE][+-]?[0-9]+)?\$/) { \$i="6.84e-7";} } print }' "${dup}/\$fname" > "\$fname"
    
                    python -u ./gp.py "\$fname" ${kary} --cadd 15 20 --gene-bed ${gene_bed} --cre-flank ${params.cre_flank}
                else
                    echo "Processing \$file from chunk folder"
                    awk -F"\\t" '{OFS=FS;for (i=7;i<=16;i++) { if(length(\$i)<1 || \$i ~ /^0+([.0]+)?([eE][+-]?[0-9]+)?\$/) { \$i="6.84e-7";} } print }' "\$file" > "\$fname"
    
                    python -u ./gp.py "\$fname" ${kary} --cadd 15 20 --gene-bed ${gene_bed} --cre-flank ${params.cre_flank}
                    
                fi
            fi
        done
    """
}
