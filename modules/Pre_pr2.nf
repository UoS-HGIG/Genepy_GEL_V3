/*
    Module: Pre_processing_2
    Purpose:
        Extract per-allele variant annotations from the VEP CSQ field and
        prepare the information required to build per-gene matrices.

    Inputs:
        - Genotype-QC-filtered VCF (f5.vcf.gz)
        - Header metadata file (header_meta)
        - GENCODE gene-region BED file (Genecode_p50_bed)
        - Pipeline templates directory, containing pre_1g_v2.sh
        - Shard number and subshard number

    Outputs:
        - Annotation columns c1-c5, plus c5a and c5b
        - Gene list (gene.lst)
        - Deduplicated VCF (f5_dedup.vcf.gz)
        - Header metadata and intermediate files used by later steps
        - CADD-missing log (cadd_missing.log)

    Key steps:
        1. Copy the annotation script `pre_1g_v2.sh` from the templates
           directory into the working directory.
        2. Prepare the metadata header and GENCODE gene-region BED file.
        3. Extract the genotype columns from the VCF.
        4. Run `pre_1g_v2.sh` to parse the VEP CSQ block and extract:
           - Locus information
           - Variant consequence
           - Gene name
           - Allele frequency
           - CADD raw score

    Notes:
        - Most of the annotation logic is implemented in
          `templates/pre_1g_v2.sh`, not in this Nextflow module.
        - This module is a key point to review when changing VEP CSQ fields,
          gnomAD frequency fields, or CADD annotation definitions.
*/
process Pre_processing_2 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  //maxForks 10
  tag "Pre_processing_2_${shard_num}_${subshard_num}"
  label "Pre_processing_2"
  

  input:
  tuple file("f5.vcf.gz"),val(shard_num),val(subshard_num) 
  path(header_meta)
  path(Genecode_p50_bed)
  path(template)
  output:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("gene.lst"),path("f5_dedup.vcf.gz"),path("header_meta"), val(shard_num),val(subshard_num) ,emit: main 
  path("f6")
  path("p1")
  path("alt")
  path("cadd_missing.log")
  //path("csq")
  //path("p1_s")
  //path("p1_m")
  //path("p1_order_1")
  //path("csq_re")
  //path("alt_re")
 // path("p1_re")
  //path("p1_1")
  //path("p1_2")
  //path("p1_order")
  //path("p1_u")
  path("c_u")
  path("c2")
  path("f61.vcf")
  path("p1.bed")
  //path("order")


  
  shell:
    """
    REAL_PATH1=\$(readlink -f ${template})
    ls \$REAL_PATH1
    cp \$REAL_PATH1/pre_pro2.sh ./pre_pro2.sh
    chmod +x ./pre_1g_v2.sh
    cat ${header_meta} > meta_CADD_head
    cat ${Genecode_p50_bed} > p50.bed
    ## bgzip -c "f5.vcf" > f5.vcf.gz
    bcftools view -h f5.vcf.gz --threads $task.cpus | grep -v "##" | cut -f 10- >p
    ./pre_pro2.sh
    """
}
