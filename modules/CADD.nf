process CADD_score {
  tag "CADD_score_${shard_num}_${subshard_num}"
  label "CADD_score"
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  maxForks 20
  input:
  tuple val(shard_num), val(subshard_num), val(chr_name), path(vcf_File), path(annotations_cadd)
 
      
 
  output:
  tuple val(shard_num),path("p1.vcf"), path("wes_${subshard_num}.raw.tsv.gz"), path("wes_${subshard_num}.raw.tsv.gz.tbi"), val(subshard_num), path(vcf_File),val(chr_name), emit: pre_proc_1
  path("${subshard_num}.p11.vcf.gz")
 
  script:
    """
    echo "CADD"
    REAL_PATH1=\$(readlink -f ${annotations_cadd})
    ln -sf \$REAL_PATH1 /opt/CADD-scripts-CADD1.6/data/annotations/GRCh38_v1.6
    
    bcftools view -G ${vcf_File} -Ov  --threads $task.cpus -o p1.vcf
    st=\$(awk '\$0 !~ /^#/ {print NR; exit}' p1.vcf)
    awk -F'\t' '$0 !~ /^#/ {print $1; exit}' p1.vcf > cadd_chr.txt
    awk -F"\t" -v OFS="\t" -v st="\$st" '
      NR < st {print; next}
      \$1 ~ /^#/ {print; next}
      { sub(/^chr/, "", \$1); print }
      ' p1.vcf  > "${subshard_num}.p11.vcf"
    bgzip -f "${subshard_num}.p11.vcf"
    tabix -p vcf "${subshard_num}.p11.vcf.gz"
    bcftools view -v indels -O v -o "${subshard_num}_indels.p11.vcf" "${subshard_num}.p11.vcf.gz"
    CADD.sh -c $task.cpus -o wes_${subshard_num}.raw.tsv.gz "${subshard_num}_indels.p11.vcf"
    tabix -p vcf wes_${subshard_num}.raw.tsv.gz
    ### CADD TSV columns: Chrom, Pos, Ref, Alt, RawScore, PHRED (col 6)
    
    """
}
