process site_qc_cadd15 {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "site_qc_cadd15_${shard_num}_${subshard_num}"
  label "site_qc_cadd15"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile)
  path(base_site_qc)
  output:
   tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile),path("combined.bed")
  script:
  
  
    """
##recommend
   echo "filter we"
   zcat "wes.tsv.gz" | awk -v OFS='\t' '
\$0 !~ /^#/ && \$6 >= 15 { print "chr"\$1, \$2-1, \$2 }' | sort -k1,1 -k2,2n > "wes_cadd15.bed"
#### tabix -s1 -b2 -e2 "wes_cadd15.bed.gz"
echo "filter WG.tsv.gz"
zcat ${plugin2} | awk -v OFS='\t' '
\$0 !~ /^#/ && \$6 >= 15 { print "chr"\$1, \$2-1, \$2 }'| sort -k1,1 -k2,2n  > "cadd15_regions.bed"
###### tabix -s1 -b2 -e2 "cadd15_regions.bed.gz"
bcftools query "${base_site_qc}/shard-${shard_num}/subshard-${subshard_num}/dragen.gel.siteqc.vcf.gz" -i '(MEDIAN_DP>=8) & (MEDIAN_GQ>=10) & (MISSINGNESS_RATE<=0.12)' -f '%CHROM\t%POS0\t%POS\n'  > "siteqc_pass_variants.bed"
echo "bcftools query done"
cat "wes_cadd15.bed" "cadd15_regions.bed" > "all_regions.bed"

echo "All region bed"

awk 'BEGIN{OFS="\t"}
NR==FNR {
    if (\$0 !~ /^#/ && \$0 != "") {
        n[\$1]++
        s[\$1,n[\$1]] = \$2
        e[\$1,n[\$1]] = \$3
    }
    next
}
{
    for (i=1; i<=n[\$1]; i++) {
        if (\$2 < e[\$1,i] && \$3 > s[\$1,i]) {
            print
            break
        }
    }
}' "siteqc_pass_variants.bed" "all_regions.bed" > "combined.bed"


    """
}