process CADD_score {
  tag "CADD_score_${shard_num}_${subshard_num}"
  label "CADD_score"
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  maxForks 20
  input:
  tuple val(shard_num), val(subshard_num), val(chr_name), path(vcf_File), path(annotations_cadd)
  // tuple val(shard_num), val(subshard_num), val (chr),  path(vcf_File), path(gnomad_joint_vcf), path(cadd_)
      
  //val cadd_param = params.cadd_
  output:
  tuple val(shard_num), path("p1.vcf"), path("wes_${subshard_num}.tsv.gz"), path("wes_${subshard_num}.tsv.gz.tbi"), val(subshard_num), path(vcf_File), emit: pre_proc_1
  path("${subshard_num}.p11.vcf")
  // path("f3b.vcf")
  script:
    """
    echo "CADD"
    REAL_PATH1=\$(readlink -f ${annotations_cadd})
    ln -sf \$REAL_PATH1 /opt/CADD-scripts-CADD1.6/data/annotations/GRCh38_v1.6
    
    bcftools view -G ${vcf_File} -Ov  --threads $task.cpus -o p1.vcf
    st=\$(awk '\$0 !~ /^#/ {print NR; exit}' p1.vcf)
    ## awk -F"\t" '\$1 ~ /^#/ || length(\$4)>1 || length(\$5)>1' p1.vcf | sed "\${st},\\\$s/chr//g" > ${subshard_num}.p11.vcf

    awk -F"\t" -v OFS="\t" -v st="\$st" '
      NR < st {print; next}
      \$1 ~ /^#/ {print; next}
      { sub(/^chr/, "", \$1); print }
      ' p1.vcf  > "${subshard_num}.p11.vcf"

    CADD.sh -c $task.cpus -o wes_${subshard_num}.tsv.gz ${subshard_num}.p11.vcf
    
    tabix -p vcf wes_${subshard_num}.tsv.gz  
    """
}
