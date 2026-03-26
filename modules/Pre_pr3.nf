process Pre_processing_3 {
  publishDir "${params.outDir}/${params.chr}/${vcf_n}", mode: "copy", overwrite: true
  //maxForks 10
  tag "Pre_processing_3_${vcf_n}"
  label "Pre_processing_3"
  input:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("gene.lst"),path("f5_dedup.vcf.gz"),path("header_meta"), val(vcf_n) , val(chrx) 
  path(template)
  output:
  path("metafiles15_*"), emit: meta_files15
  path("metafiles20_*"), emit: meta_files20
  path("metafilesALL_*"), emit: meta_filesALL
  //tuple path("metafilesALL"),path("metafiles15"),path("metafiles20"), emit: folders
  //path("*")
  //file("c6")
  //file("meta_CADDALL.txt")
  shell:
    """
    echo "Processing 3"
    echo "vcf_n: ${vcf_n}"

    REAL_PATH1=\$(readlink -f ${template})
    ##ls \$REAL_PATH1
    cp \$REAL_PATH1/pre_2.sh ./pre_2.sh
    chmod +x ./pre_2.sh
    VCF_NAME=\$(basename ${vcf_n})
    region=\$(echo \$VCF_NAME | awk -F'[_|.]' '{print \$5"_"\$6}')
    ##mkdir -p metafilesALL metafiles15 metafiles20
    ##touch metafilesALL/ALL.txt
    ##touch metafiles15/15.txt
    ##touch metafiles20/20.txt
    ./pre_2.sh \$region
    """
}
