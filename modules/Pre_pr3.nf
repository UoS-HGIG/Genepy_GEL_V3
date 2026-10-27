/*
    Module: Pre_processing_3
    Purpose:
        Convert the extracted variant annotations and genotypes into per-gene
        `.meta` input matrices for GenePy scoring.

    Inputs:
        - Annotation outputs from Pre_processing_2:
          c1, c2, c3, c4, c5, c5a, c5b, gene.lst, f5_dedup.vcf.gz,
          and header metadata
        - Templates directory, containing pre_2.sh
        - Shard number and subshard number

    Outputs:
        - Per-gene `.meta` matrices for variants with CADD PHRED >= 15
          (metafiles15_*)
        - Per-gene `.meta` matrices for variants with CADD PHRED >= 20
          (metafiles20_*)
        - Summary matrices meta_CADD15.txt and meta_CADD20.txt

    Key steps:
        1. Copy `pre_2.sh` from the templates directory into the task
           working directory.
        2. Run `pre_2.sh` with the shard and subshard identifiers.
        3. Generate separate per-gene matrices at two CADD thresholds:
           15 and 20.

    Notes:
        - The `.meta` files contain one row per variant, with annotation columns
          followed by genotype columns.
        - These matrices are the direct inputs to GenePy scoring.
        - The core matrix-generation logic is implemented in
          `templates/pre_2.sh`.
*/
process Pre_processing_3 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  //maxForks 10
  tag "Pre_processing_3_${shard_num}_${subshard_num}"
  label "Pre_processing_3"
  input:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("gene.lst"),path("f5_dedup.vcf.gz"),path("header_meta"), val(shard_num),val(subshard_num)
  path(template)
  output:
  path("metafiles15_*"), emit: meta_files15
  path("metafiles20_*"), emit: meta_files20
  path("meta_CADD15.txt")
  path("meta_CADD20.txt")
  shell:
    """
    echo "Processing 3"
    echo "vcf_n: ${shard_num}_${subshard_num}"

    REAL_PATH1=\$(readlink -f ${template})
    cp \$REAL_PATH1/pre_2.sh ./pre_pro3.sh
    chmod +x ./pre_2.sh
    ./pre_pro3.sh "${shard_num}" "${subshard_num}"
    """
}
