process Pre_processing_1 {
  publishDir "${params.outDir}/${params.chr}/${vcf_n}", mode: "copy", overwrite: true
  // maxForks 10
  tag "Pre_processing_1_${vcf_n}"
  label "Pre_processing_1"
  
  input:
  tuple path(x), val(vcf_n), path(vcfFile), val(chrx),path("input.vcf.gz")
  path(ethnicity)
  path(xgen_bed)
  output:
  tuple path("f5.vcf.gz"), val(vcf_n), val(chrx), emit:main
  path("*.vcf.gz")
  path("*.vcf")
  
  
  shell:
    """
    gunzip -c "input.vcf.gz" | grep -v '##'|cut -f 9-> p2
    grep -v '##' ${x} > p1
    grep '##' ${x} > f31.vcf
    paste p1 p2 >> f31.vcf
    rm -r p1 p2
    awk -F"\t" '\$7~/PASS/ || \$1~/#/' f31.vcf > f3.vcf
    bcftools view -h  "input.vcf.gz" --threads $task.cpus | grep '^##FORMAT=' > format.txt
    sed -i '1 r format.txt' f3.vcf
    #####
    ##bcftools +fill-tags f3.vcf --threads $task.cpus -- -t 'FORMAT/AB:1=float((AD[:1]) / (DP))' | bgzip -c > f3.vcf.gz
    ##rm f3.vcf
    ##tabix -p vcf f3.vcf.gz
    ##bcftools filter -S . --include '(FORMAT/DP>=8 & FORMAT/AB>=0.15) |FORMAT/GT="0/0" | FORMAT/GT="0"'  --threads $task.cpus -Oz -o f3b.vcf.gz f3.vcf.gz
    ### bcftools filter -S . --include 'FORMAT/DP>=8 & FORMAT/AB>=0.15 |FORMAT/GT="0/0"'  --threads $task.cpus -Oz -o f3b.vcf.gz f3.vcf.gz
    ##tabix -p vcf f3b.vcf.gz

    cat ${ethnicity} > ethnicity.txt
    bcftools +fill-tags f3.vcf --threads $task.cpus -- -S ethnicity.txt -t 'HWE,F_MISSING' | bcftools view -e '(CHROM=="chrY" & INFO/F_MISSING>=0.56 & INFO/HWE_1>(0.05/15922704))' --threads $task.cpus -Ov -o f4.vcf 
    bcftools view -i 'INFO/F_MISSING<0.12 & INFO/HWE_1>(0.05/15922704)' --threads $task.cpus -Oz -o f5.vcf.gz f4.vcf

    tabix -p vcf f5.vcf.gz
    
    ##bcftools view f4.vcf.gz -R ${xgen_bed} --threads $task.cpus -Oz -o f5.vcf.gz
    
    
    """
}
