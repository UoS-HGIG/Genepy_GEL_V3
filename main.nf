#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.shard_path
]
 
include { CADD_score } from "./modules/CADD"  
include { site_qc_cadd15 } from "./modules/site_qc_cadd15"
include { VEP_score } from "./modules/VEP"  
include { Pre_processing_1 } from "./modules/Pre_pr1"  
include { Pre_processing_2 } from "./modules/Pre_pr2"  
include { Pre_processing_3 } from "./modules/Pre_pr3"  
include { Reatt_Genes } from "./modules/Gene_reattach"
include { Genepy_score } from "./modules/Genepy"



// Define workflow
workflow {

    println """\
         G E N E P Y           P I P E L I N E
          ===================================
     G E N O M I C --------------- M E D I C I N E 
                         UoS
                     Iman Nazari
          ===================================
          Samples         : ${params.shard_path}
         """.stripIndent()
     
        def shard_dir_name = file(params.shard_path).name
        
        def target_chr = params.chr
        def shard_path_pattern = "${params.shard_path}/shard-*/subshard-*/dragen.vcf.gz"
       println "START :)"
      


        def shard_map = [:]

        file(params.multiallelic_shards_bed)
    .eachLine { line ->
        if( !line?.trim() ) return

        def cols = line.split('\\t| +')
        def chr_name = cols[0]
        def shard    = cols[4].toString()
        def subshard = cols[5].toString()

        if( chr_name == target_chr ) {
            shard_map["${shard}_${subshard}"] = chr_name
        }
    }

chrx = Channel.fromPath(shard_path_pattern, checkIfExists: true)
.filter { vcf_file ->
    vcf_file.parent.parent.name == "shard-93" &&
    vcf_file.parent.name == "subshard-18"
}
    .map { vcf_file ->
        def shard_num       = vcf_file.parent.parent.name.replace('shard-', '')
        def subshard_number = vcf_file.parent.name.replace('subshard-', '')
        def chr_name        = shard_map["${shard_num}_${subshard_number}"]

        if( chr_name ) {
            tuple(shard_num, subshard_number, chr_name, vcf_file, file(params.annotations_cadd))
        } else {
            null
        }
    }
    .filter { it != null }
  
      CADD_score(chrx)
      site_qc_cadd15(CADD_score.out.pre_proc_1,params.base_site_qc,params.plugin3)
//      qc_split = site_qc_cadd15.out.site_qc_out.branch { shard_num, p1, wes, tbi, subshard_num, vcfFile, combined_bed ->
//        has_variants: combined_bed.size() > 0
//        empty: combined_bed.size() == 0
//    }
      VEP_score(site_qc_cadd15.out.site_qc_out,params.homos_vep,params.vep_plugins,params.plugin1,params.plugin2,params.genomad_indx1,params.genomad_indx2,params.base_site_qc)
      Pre_processing_1(VEP_score.out.vep_out)
      Pre_processing_2(Pre_processing_1.out.main,params.header_meta,params.gene_code_bed,params.templates)
      Pre_processing_3(Pre_processing_2.out.main,params.templates)     
      def meta15 = Pre_processing_3.out.meta_files15.collect().map { genes_list -> ["15",params.chr, genes_list] }
      def meta20 = Pre_processing_3.out.meta_files20.collect().map { genes_list -> ["20",params.chr, genes_list] }

      x_combo= meta15.concat(meta20)
      Reatt_Genes(x_combo)
   

def flatMetas = Reatt_Genes.out.path.flatten()
def flatDups  = Reatt_Genes.out.dup.flatten()


def metas = flatMetas.map { p ->
    def fullKey = p.toString().tokenize('/').find { it.startsWith('metafiles') }
    def baseKey = fullKey.replaceAll(/(_\d+)+$/, '')  
    tuple(baseKey, p)                                
}


def dups = flatDups.map { d ->
    def fullKey = d.toString().tokenize('/').find { it.startsWith('dup') }
    def baseKey = fullKey?.replace('dup', 'metafiles')
    tuple(baseKey, d)                                
}
dups.collect()

def met_ = metas
    .combine(dups)                       // produce all pairs
    .filter { m -> m[0] == m[2] }     // keep only matches on baseKey
    .map { m ->                       // 
        def key         = m[0]
        def folder_path = m[1]
        def dup_path    = m[3]

        def cadd_score = (key == 'metafiles20') ? '20' :
                         (key == 'metafiles15') ? '15' : '15'

        tuple(folder_path, params.chr, cadd_score, params.genepy_py,params.kary, dup_path)
    }
    .view()
     Genepy_score(met_)
}
workflow.onComplete {
   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}




