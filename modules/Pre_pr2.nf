process Pre_processing_2 {
  publishDir "${params.outDir}/${params.chr}/${vcf_n}", mode: "copy", overwrite: true
  //maxForks 10
  tag "Pre_processing_2_${vcf_n}"
  label "Pre_processing_2"
  

  input:
  tuple file("f5.vcf.gz"), val(vcf_n), val(chrx) 
  path(header_meta)
  path(IBD_gwas_bed)
  path(Genecode_p50_bed)
  path(template)
  output:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("gene.lst"),path("f5_dedup.vcf.gz"),path("header_meta"), val(vcf_n) , val(chrx) ,emit: main 
  
  shell:
    """
    REAL_PATH1=\$(readlink -f ${template})
    ls \$REAL_PATH1
    cp \$REAL_PATH1/pre_1.sh ./pre_1.sh
    chmod +x ./pre_1.sh
    cat ${header_meta} > meta_CADD_head
    cat ${IBD_gwas_bed} > IBD.bed
    cat ${Genecode_p50_bed} > p50.bed
    ## bgzip -c "f5.vcf" > f5.vcf.gz
    bcftools view -h f5.vcf.gz --threads $task.cpus | grep -v "##" | cut -f 10- >p
    ./pre_1.sh
    """
}
