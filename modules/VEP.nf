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
  
  output:
   tuple path("${subshard_num}.p1.vep.vcf.gz"), file(vcfFile), val(shard_num) ,emit: vep_out
   path("${subshard_num}.full_cadd15.vcf.gz")
  script:
  
  
    """
##recommend
   echo "retain cadd>=15"
   zcat "wes.tsv.gz" | awk -v OFS='\t' '
NR==1 { print; next }
\$0 !~ /^#/ && \$6 >= 15 { print }
' | bgzip -c > "wes_cadd15.tsv.gz"


    vep -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADe_AF,CADD_PHRED" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
##vep -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --max_af  --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,MAX_AF,MAX_AF_POPS,CADD_RAW,CADD_PHRED" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
  echo "2" 
  bcftools +split-vep \
      -c "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADe_AF,CADD_PHRED:Float" \
      -s worst \
      -i 'CADD_PHRED >= 15 || CADD_PHRED = "."' \
      "${subshard_num}.p1.vep.vcf" \
      -O z -o "${subshard_num}.CADD15.p1.vep.vcf.gz"

  tabix -p vcf "${subshard_num}.CADD15.p1.vep.vcf.gz"
  ##bgzip -c "${subshard_num}.p1.vep.vcf" > "${subshard_num}.p1.vep.vcf.gz"
  ##tabix -p vcf "${subshard_num}.p1.vep.vcf.gz"
  tabix -p vcf "${vcfFile}"
echo "3"
    bcftools isec -n=2 -w1 -c both -Oz  -o "${subshard_num}.full_cadd15.vcf.gz" ${vcfFile}  "${subshard_num}.CADD15.p1.vep.vcf.gz"

    tabix -p vcf "${subshard_num}.full_cadd15.vcf.gz"
    """
}
