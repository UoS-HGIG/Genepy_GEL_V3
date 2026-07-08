process site_qc_cadd15 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "site_qc_cadd15_${shard_num}_${subshard_num}"
  label "site_qc_cadd15"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile),val(chr_name)
  path(base_site_qc)
  path(plugin3)
  output:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile),path("combined.bed"), emit: site_qc_out
  path("has_variants.txt"), emit: has_variants
  path("all_regions.bed")
  path("siteqc_pass_variants.bed")
  path("wes_cadd15.bed")
  script:
  
  
    """
##recommend
   echo "filter we"
   zcat "wes.tsv.gz" | awk -v OFS='\\t' '
\$0 !~ /^#/ && \$6 >= 15 { print "chr"\$1, \$2-1, \$2 }' | sort -k1,1 -k2,2n > "wes_cadd15.bed"
#### tabix -s1 -b2 -e2 "wes_cadd15.bed.gz"
echo "filter WG.tsv.gz"
######zcat ${plugin3} | awk -v OFS='\t' '
######\$0 !~ /^#/ && \$6 >= 15 { print "chr"\$1, \$2-1, \$2 }'| sort -k1,1 -k2,2n  > "cadd15_regions.bed"
###### tabix -s1 -b2 -e2 "cadd15_regions.bed.gz"
if [[ "${chr_name}" == "chrX" || "${chr_name}" == "X" ]]; then
      echo "chrX detected: using chrX-specific siteqc filter"
      bcftools query "${base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz" \
        -i '((MEDIAN_DP_XX>=8) && (MEDIAN_GQ_XX>=10) && (MISSINGNESS_RATE_XX<=0.12))||((MEDIAN_DP_XY>=8) && (MEDIAN_GQ_XY>=10) && (MISSINGNESS_RATE_XY<=0.12))' \
        -f '%CHROM\\t%POS0\\t%POS\n' > "siteqc_pass_variants.bed"
  else
      echo "autosome/other chromosome detected: using default siteqc filter"
      bcftools query "${base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz" \
        -i '(MEDIAN_DP>=8) && (MEDIAN_GQ>=10) && (MISSINGNESS_RATE<=0.12)' \
        -f '%CHROM\\t%POS0\\t%POS\n' > "siteqc_pass_variants.bed"
  fi
echo "bcftools query done"
head "siteqc_pass_variants.bed"
head "wes_cadd15.bed"
head ${plugin3}

cat "wes_cadd15.bed" ${plugin3} > "all_regions.bed"
head "all_regions.bed"
echo "All region bed"
awk 'BEGIN{OFS="\\t"} \$1 !~ /^#/ && NF>=3 {print \$1, \$2, \$3}' "all_regions.bed" > "all_regions.clean.bed"
head "all_regions.clean.bed"
bedtools intersect -u -a "all_regions.clean.bed"  -b "siteqc_pass_variants.bed"  > "combined.bed"
if [ ! -s "combined.bed" ]; then
    echo 0 > has_variants.txt
elif [ "$(grep -vc '^#' "combined.bed" || true)" -eq 0 ]; then
    echo 0 > has_variants.txt
else
    echo 1 > has_variants.txt
fi

    """
}
