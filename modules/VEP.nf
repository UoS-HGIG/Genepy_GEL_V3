process VEP_score {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "VEP_score_${shard_num}_${subshard_num}"
  label "VEP_score"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile),path(gnomad_joint_vcf)
  path(homos_vep)
  path(vep_plugins)
  path(plugin1)
  path(plugin2)
  path(genomad_indx1)
  path(genomad_indx2)
  
  output:
   tuple path("${subshard_num}.p1.vep.vcf"), file(vcfFile), val(shard_num) ,emit: vep_out
  
  script:
  
    """

    vep -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg --max_af --custom "${gnomad_joint_vcf}",gnomad_joint,vcf,exact,0,AF_grpmax_joint --fields Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,MAX_AF,MAX_AF_POPS,AF_grpmax_joint,CADD_RAW,CADD_PHRED" -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
    """
}
