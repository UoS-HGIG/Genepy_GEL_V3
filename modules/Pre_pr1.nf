process Pre_processing_1 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 10
  tag "Pre_processing_1_${shard_num}_${subshard_num}"
  label "Pre_processing_1"
  
  input:
  tuple path("full_cadd15.vcf.gz"), val(shard_num),val(subshard_num)
  ////path("vep_out")
  output:
  tuple path("${shard_num}_${subshard_num}_f3.vcf.gz"), val(shard_num),val(subshard_num), emit:main
  path("f3_1.vcf.gz")
 
  script:
    """
echo "step1"
tabix -p vcf "full_cadd15.vcf.gz"
bcftools +fill-tags "full_cadd15.vcf.gz" -- -t 'FORMAT/AB:1=float((FORMAT/LAD[:1]) / (FORMAT/DP))' | bgzip -c > f3_1.vcf.gz
tabix -p vcf f3_1.vcf.gz
####
echo "step2"
bcftools filter -S . --include '(FORMAT/FT="PASS") && ((FORMAT/DP>=8 & FORMAT/AB>=0.15) || (FORMAT/GT="0/0") || (FORMAT/GT="0"))' -Oz -o "${shard_num}_${subshard_num}_f3.vcf.gz" f3_1.vcf.gz
   
    
    
    """
}
