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
   tuple path("${subshard_num}.p1.vep.vcf"), file(vcfFile), val(shard_num) ,emit: vep_out
   path("${subshard_num}.full_cadd15.vcf.gz")
  script:
  
  
    """
##recommend
    bgzip -c "p1.vcf" > "p1.vcf.gz"
    tabix -p vcf "p1.vcf.gz"
    zcat "wes.tsv.gz" | awk 'NR>2 && \$1 !~ /^#/ {print \$1"\t"\$2}' > cadd15_regions.txt
    bcftools view -R cadd15_regions.txt -Ov "p1.vcf.gz" > p1.cadd15.vcf
    vep -i p1.cadd15.vcf --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADe_AF,CADD_RAW" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
##vep -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --max_af  --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,MAX_AF,MAX_AF_POPS,CADD_RAW,CADD_PHRED" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
    bgzip -c "${subshard_num}.p1.vep.vcf" > "${subshard_num}.p1.vep.vcf.gz"
    tabix -p vcf "${subshard_num}.p1.vep.vcf.gz"
    tabix -p vcf "${vcfFile}"
    bcftools isec -n=2 -w1 -Oz -o "${subshard_num}.full_cadd15.vcf.gz" ${vcfFile} "${subshard_num}.p1.vep.vcf.gz"

    tabix -p vcf "${subshard_num}.full_cadd15.vcf.gz"
    """
}
