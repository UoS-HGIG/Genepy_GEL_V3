/*
    Module: Pre_processing_1
    Purpose:
        Perform sample-level genotype QC on the VEP-annotated,
        CADD >= 15 genotype VCF.

    Inputs:
        - Genotype-bearing VCF with VEP CSQ annotations (full_cadd15.vcf.gz)
        - Shard number and subshard number

    Outputs:
        - Genotype-QC-filtered VCF (<shard>_<subshard>_f3.vcf.gz)
        - Intermediate VCF with allele balance added (f3_1.vcf.gz)

    Key steps:
        1. Calculate FORMAT/AB (allele balance) as:
           LAD first value divided by DP.
        2. Set genotypes to missing when they fail the following criteria:
           - FORMAT/FT must be PASS.
           - Either DP >= 8 and AB >= 0.15, or the genotype is homozygous
             reference (0/0 or 0).
        3. Compress and index the filtered VCF.

    Notes:
        - Homozygous-reference genotypes are retained even if they do not meet
          the depth and allele-balance thresholds.
        - This module filters genotypes within variants; it does not remove
          variant sites.
*/
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
