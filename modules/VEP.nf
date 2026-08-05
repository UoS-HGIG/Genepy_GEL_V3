process VEP_score {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "VEP_score_${shard_num}_${subshard_num}"
  label "VEP_score"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile), path("combined1.bed")
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
   path("combined.bed")
  script:
  
  
    """


bgzip -c "p1.vcf" > "p1.vcf.gz"
tabix -p vcf "p1.vcf.gz"
echo "intersect p1 & WG"
sort "combined1.bed" | uniq > "combined.bed"
bcftools view -R "combined.bed" -O z  --threads $task.cpus -o "filtered_cadd15.vcf.gz" "p1.vcf.gz"
tabix -p vcf "filtered_cadd15.vcf.gz"

echo "VEP start"
    vep -i "filtered_cadd15.vcf.gz" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${plugin1},${plugin2},"wes.tsv.gz" --af_gnomade --af_gnomadg -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
echo "VEP done"
 

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
