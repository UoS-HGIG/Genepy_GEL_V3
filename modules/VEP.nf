/*
    Module: VEP_score
    Purpose:
        Annotate site-QC-passing, CADD >= 15 variants with VEP consequence,
        gnomAD allele-frequency, and CADD annotations.

    Inputs:
        - Filtered CADD >= 15 VCF from the site-QC module
        - Original genotype-bearing subshard VCF
        - VEP cache directory
        - VEP plugin directory
        - Precomputed CADD SNV scores (plugin1)
        - Precomputed CADD indel scores (plugin2)
        - gnomAD index files required by the VEP CADD plugin

    Outputs:
        - VEP-annotated sites-only VCF (<subshard>.p1.vep.vcf.gz)
        - Genotype-bearing VCF with CSQ annotations added
          (<subshard>.full_cadd15.vcf.gz)

    Key steps:
        1. Run VEP offline against GRCh38 using the configured VEP cache.
        2. Use the CADD VEP plugin, supplying precomputed SNV and indel CADD
           scores plus the CADD indel output from the earlier CADD module.
        3. Add gnomAD exomes and genomes allele-frequency annotations using
           `--af_gnomade` and `--af_gnomadg`.
        4. Intersect the VEP-annotated sites-only VCF with the original
           genotype-bearing VCF.
        5. Transfer the VEP CSQ annotation block onto the genotype-bearing VCF
           using `bcftools annotate -c '+CSQ'`.

    Notes:
        - The output CSQ field contains the annotations parsed in the following
          preprocessing modules.
        - CADD scores for SNVs come from the precomputed whole-genome SNV file;
          CADD scores for indels come from the CADD module.
*/
process VEP_score {
  publishDir "${params.outDir}/${shard_num}/${subshard_num}", mode: "copy", overwrite: true
  // maxForks 20
  tag "VEP_score_${shard_num}_${subshard_num}"
  label "VEP_score"
  
  input:
  tuple val(shard_num), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(subshard_num) , path(vcfFile), path("filtered_cadd15.vcf.gz")
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
  script:
  
  
    """
  tabix -p vcf "filtered_cadd15.vcf.gz"
echo "VEP start"
  vep -i "filtered_cadd15.vcf.gz" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --allele_number --plugin CADD,"${plugin2}","wes.tsv.gz" --af_gnomade --af_gnomadg -o "${subshard_num}.p1.vep.vcf" --dir_cache ${homos_vep}  --dir_plugins ${vep_plugins}
echo "VEP done"
   bgzip -c "${subshard_num}.p1.vep.vcf" > "${subshard_num}.p1.vep.vcf.gz"
   tabix -p vcf "${subshard_num}.p1.vep.vcf.gz"
   tabix -p vcf "${vcfFile}"
 echo "adding vep  annotation output to the original vcf file"
    bcftools isec -n=2 -w1 -c both -Oz  -o "${subshard_num}.isec_cadd15.vcf.gz" ${vcfFile}  "${subshard_num}.p1.vep.vcf.gz"
    tabix -p vcf "${subshard_num}.isec_cadd15.vcf.gz"
    bcftools annotate -a   "${subshard_num}.p1.vep.vcf.gz" -c '+CSQ'  -Oz -o "${subshard_num}.full_cadd15.vcf.gz" "${subshard_num}.isec_cadd15.vcf.gz"
    bcftools index -f "${subshard_num}.full_cadd15.vcf.gz"
    """
}
