process Pre_processing_1 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 10
  tag "Pre_processing_1_${shard_num}_${subshard_num}"
  label "Pre_processing_1"
  
  input:
  tuple path("full_cadd15.vcf.gz"), val(shard_num),val(subshard_num)
  path("vep_out")
  output:
  tuple path("${shard_num}_{subshard_num}_f3.vcf.gz"), val(shard_num),val(subshard_num), emit:main
  path("f3_1.vcf.gz")
  
  script:
    """
#### add cadd filteration >=15  and make sure datframe is correct ( do we need normalization or nor)
######GC comment####
####step 1: variant/site based filtration
##with CADD15 filtration and filter for varaints with QUAL="PASS";  the following is site-qc filtration for autosome variant
##the input file 'dragen.gel.siteqc.vcf.gz' for each vcf can be found in the following file
  ##s3://512426816668-gel-data-resources/dragen3.7.8/AggV3_resources/manifests/site_qc/2026-01-06/siteqc_shards.bed

##bcftools query filesystems/dragen.gel.siteqc.vcf.gz -i '(MEDIAN_DP>=8) & (MEDIAN_GQ>=10) & (MISSINGNESS_RATE<=0.12)' -f '%CHROM:%POS:%REF:%ALT\n' > siteqc_pass_variants.tsv
##bcftools view -i 'ID=@siteqc_pass_variants.tsv' filesystems/dragen.vcf.gz -Oz -o pass_variants_filtered.vcf.gz

##after this is genotype-based quality control
##bcftools filter -S . --include '(FORMAT/DP>=8 & FORMAT/AB>=0.15) |FORMAT/GT="0/0" | FORMAT/GT="0"'  --threads $task.cpus -Oz -o pass_variants_filtered.vcf.gz f3.vcf.gz
########End of QC and variant filtration for pre-processing########
echo "step1"
bcftools query "${params.base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz" -i '(MEDIAN_DP>=8) & (MEDIAN_GQ>=10) & (MISSINGNESS_RATE<=0.12)' -f '%CHROM:%POS:%REF:%ALT\n' > siteqc_pass_variants.tsv
echo "step2"
cat > tmp.vcf <<'EOF'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
EOF
echo "step3"
zcat "vep_out" | grep "#CHROM" >> tmp.vcf

awk -F':' 'BEGIN{OFS="\t"}
NF>=4 {
    print \$1, \$2, ".", \$3, \$4, ".", "PASS", "."
}' siteqc_pass_variants.tsv >> tmp.vcf

bgzip -f tmp.vcf
tabix -f -p vcf tmp.vcf.gz
echo "step4"
bcftools isec -n=2 -w1 -Oz -o siteqc_pass_variants_filtered.vcf.gz "full_cadd15.vcf.gz" tmp.vcf.gz
bcftools +fill-tags siteqc_pass_variants_filtered.vcf.gz -- -t 'FORMAT/AB:1=float((FORMAT/AD[:1]) / (FORMAT/DP))' | bgzip -c > f3_1.vcf.gz
tabix -p vcf f3_1.vcf.gz
####
echo "step5"
bcftools filter -S . --include '(FORMAT/FT="PASS") && ((FORMAT/DP>=8 & FORMAT/AB>=0.15) | FORMAT/GT="0/0" | FORMAT/GT="0")' -Oz -o ${shard_num}_{subshard_num}_f3.vcf.gz f3_1.vcf.gz
   
    
    
    """
}
