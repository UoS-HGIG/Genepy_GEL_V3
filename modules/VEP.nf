process VEP_score {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "VEP_score_${shard_num}_${subshard_num}"
  label "VEP_score"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile)
  path(homos_vep)
  path(vep_plugins)
  path(plugin1)
  path(plugin2)
  path(genomad_indx1)
  path(genomad_indx2)
  path(base_site_qc)
  output:
   tuple path("${subshard_num}.full_cadd15.vcf.gz"), val(shard_num),val(subshard_num) ,emit: vep_out
   path("${subshard_num}.p1.vep.vcf.gz"),emit: vep_out2
   path(combined.bed)
   path("filtered_cadd15.vcf.gz")
   path("siteqc_pass_variants.bed")
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
cat "wes_cadd15.bed" "cadd15_regions.bed" > all_regions.bed

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
}' all_regions.bed "siteqc_pass_variants.bed" > combined.bed
lcout=\$(wc -l < combined.bed)
echo "combined lines $lcout"
bgzip -c "p1.vcf" > "p1.vcf.gz"
tabix -p vcf "p1.vcf.gz"
echo "intersect p1 & WG"
bcftools view -R combined.bed -O z  --threads $task.cpus -o "filtered_cadd15.vcf.gz" "p1.vcf.gz"
tabix -p vcf "filtered_cadd15.vcf.gz"

echo "VEP"
    vep -i "filtered_cadd15.vcf.gz" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
##vep -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --max_af  --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,MAX_AF,MAX_AF_POPS,CADD_RAW,CADD_PHRED" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
  echo "2" 
 ## bcftools +split-vep \
 ##     -c "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADe_AF,CADD_PHRED:Float" \
 ##     -s worst \
 ##     -i 'CADD_PHRED >= 15 || CADD_PHRED = "."' \
 ##     "${subshard_num}.p1.vep.vcf" \
 ##     -O z -o "${subshard_num}.CADD15.p1.vep.vcf.gz"

  ##tabix -p vcf "${subshard_num}.p1.vep.vcf"
  bgzip -c "${subshard_num}.p1.vep.vcf" > "${subshard_num}.p1.vep.vcf.gz"
  
  tabix -p vcf "${subshard_num}.p1.vep.vcf.gz"
  tabix -p vcf "${vcfFile}"
echo "3"
    bcftools isec -n=2 -w1 -c both -Oz  -o "${subshard_num}.isec_cadd15.vcf.gz" ${vcfFile}  "${subshard_num}.p1.vep.vcf.gz"
    tabix -p vcf "${subshard_num}.isec_cadd15.vcf.gz"
    bcftools annotate -a   "${subshard_num}.p1.vep.vcf.gz" -c '+CSQ'  -Oz -o "${subshard_num}.full_cadd15.vcf.gz" "${subshard_num}.isec_cadd15.vcf.gz"
    bcftools index -f "${subshard_num}.full_cadd15.vcf.gz"
    """
}
