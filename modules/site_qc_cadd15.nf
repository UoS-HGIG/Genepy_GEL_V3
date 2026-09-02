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
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile),path("filtered_cadd15.vcf.gz"), emit: site_qc_out
  path("all_regions.clean.bed")
  path("siteqc_pass.vcf.gz")
  path("wes_cadd15.bed")
  path("siteqc_pass.multi.vcf.gz")
  path("filtered_cadd15a.vcf.gz")
  path("siteqc_pass_multi_sites.tsv")
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
      bcftools view \
      -i '((MEDIAN_DP_XX>=8) && (MEDIAN_GQ_XX>=10) && (MISSINGNESS_RATE_XX<=0.12))||((MEDIAN_DP_XY>=8) && (MEDIAN_GQ_XY>=10) && (MISSINGNESS_RATE_XY<=0.12))' \
      -O z --threads $task.cpus \
      -o "siteqc_pass.vcf.gz" \
      "${base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz"

      tabix -p vcf "siteqc_pass.vcf.gz"
  else
      echo "autosome/other chromosome detected: using default siteqc filter"
      bcftools view \
      -i '(MEDIAN_DP>=8) && (MEDIAN_GQ>=10) && (MISSINGNESS_RATE<=0.12)' \
      -O z --threads $task.cpus \
      -o "siteqc_pass.vcf.gz" \
      "${base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz"

      tabix -p vcf "siteqc_pass.vcf.gz"
  fi
echo "bcftools filteration step has done"
zcat "siteqc_pass.vcf.gz" | tail

##convert siteqc vcf to multiallelic vcf
bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SOURCE_RECORD\\n' "siteqc_pass.vcf.gz" \\
| awk -v OFS='\\t' '{
    if (\$5 != "." && \$5 != "") { split(\$5, a, "|"); print a[1], a[2], a[3], a[4] }
    else                        { print \$1, \$2, \$3, \$4 }
  }' \\
| sort -u | sort -k1,1 -k2,2n > "siteqc_pass_multi_sites.tsv"

# the two split rows of a multi-allelic locus carry the same SOURCE_RECORD,
# so "sort -u" collapses them back to the single record held in p1
bgzip -c "p1.vcf" > "p1.vcf.gz"
tabix -p vcf "p1.vcf.gz"
bcftools view -h "p1.vcf.gz" | grep -E '^##(fileformat|contig|reference)' > "siteqc_pass.multi.vcf"
printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' >> "siteqc_pass.multi.vcf"
awk -v OFS='\\t' '{print \$1, \$2, ".", \$3, \$4, ".", ".", "."}' "siteqc_pass_multi_sites.tsv" >> "siteqc_pass.multi.vcf"

bgzip -f "siteqc_pass.multi.vcf"
tabix -p vcf "siteqc_pass.multi.vcf.gz"
echo "site QC file restored to multi-allelic representation"
##head "wes_cadd15.bed"
head ${plugin3}
bgzip -c "p1.vcf" > "p1.vcf.gz"
tabix -p vcf "p1.vcf.gz"

cat "wes_cadd15.bed" ${plugin3} > "all_regions.bed"
head "all_regions.bed"
echo "All region bed"
awk 'BEGIN{OFS="\\t"} \$1 !~ /^#/ && NF>=3 {print \$1, \$2, \$3}' "all_regions.bed" > "all_regions.clean.bed"
head "all_regions.clean.bed"

bcftools view -R "all_regions.clean.bed" -O z  --threads $task.cpus -o "filtered_cadd15a.vcf.gz" "p1.vcf.gz"
tabix -p vcf "filtered_cadd15a.vcf.gz"
bcftools isec -c none -n=2 -w1 -O z --threads $task.cpus -o "filtered_cadd15.vcf.gz" "filtered_cadd15a.vcf.gz" "siteqc_pass.multi.vcf.gz"
tabix -p vcf "filtered_cadd15.vcf.gz"

    """
}
