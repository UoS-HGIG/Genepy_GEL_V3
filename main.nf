#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.vcf
]
 
include { CADD_score } from "./modules/CADD"  
include { VEP_score } from "./modules/VEP"  
////include { Pre_processing_1 } from "./modules/Pre_pr1"  
////include { Pre_processing_2 } from "./modules/Pre_pr2"  
////include { Pre_processing_3 } from "./modules/Pre_pr3"  
////include { Reatt_Genes } from "./modules/Gene_reattach"
////include { Genepy_score } from "./modules/Genepy"



// Define workflow
workflow {

    println """\
         G E N E P Y           P I P E L I N E
          ===================================
     G E N O M I C --------------- M E D I C I N E 
                         UoS
                     Iman Nazari
          ===================================
         Samples         : ${params.vcf}
         params.cadd     : ${params.annotations_cadd}
         """.stripIndent()
     
        
        def shard_number = params.shard_number
        def shard_path_pattern = "${params.vcf}/shard-${shard_number}/subshard-*/dragen.vcf.gz"
       println "START :)"
       println "Shard list: $shard_number"


        def shard_map = [:]

        file(params.multiallelic_shards_bed)
            .eachLine { line ->
                if( !line?.trim() ) return

                def cols = line.split('\\t| +')
                def chr_name = cols[0]
                def shard    = cols[4].toString()
                def subshard = cols[5].toString()

                shard_map["${shard}_${subshard}"] = chr_name
            }   

        chrx = Channel.fromPath(shard_path_pattern, checkIfExists: true)
        .map { vcf_file->
            def shard_num       = params.shard_number.toString()
            def subshard_number = f.parent.name.replace('subshard-', '')
            def vcf_n = "subshard-${subshard_number}"
            def chr_name        = shard_map["${shard_num}_${subshard_number}"]
            if( !chr_name ) {
                throw new IllegalArgumentException("No chromosome mapping found for shard=${shard_num}, subshard=${subshard_num}")
            }

            def gnomad_joint_vcf = "${params.gnomad_joint_dir}/${chr_name}.joint.vcf.gz"
            def vcf_n = "Shard${shard_num}_Subshard_${subshard_num}"

            tuple(shard_num, subshard_num, chr_name, vcf_file, file(gnomad_joint_vcf),file(params.annotations_cadd))
        }
  //          tuple(shard_num, vcf_n, f, file(params.annotations_cadd))
  //      }
      CADD_score(chrx)
      VEP_score(CADD_score.out.pre_proc_1,params.homos_vep,params.vep_plugins,params.plugin1,params.plugin2,params.genomad_indx1,params.genomad_indx2)
////      Pre_processing_1(VEP_score.out,params.ethnicity,params.xgen_bed)
////      Pre_processing_2(Pre_processing_1.out.main,params.header_meta,params.IBD_gwas_bed,params.Genecode_p50_bed,params.templates)
////      Pre_processing_3(Pre_processing_2.out.main,params.templates)     
////     // def meta15 = Pre_processing_3.out.meta_files15.collect().map { genes_list -> ["15",chromosomeList, genes_list] }
////      def meta15 = Pre_processing_3.out.meta_files15.collect().map {genes_list -> ["15",chromosomeList, genes_list] }
////
////      def meta20 = Pre_processing_3.out.meta_files20.collect().map { genes_list -> ["20",chromosomeList, genes_list] }
////      def metaALL = Pre_processing_3.out.meta_filesALL.collect().map { genes_list -> ["ALL",chromosomeList, genes_list] }
////      x_combo= meta15.concat(meta20).concat(metaALL)
////      Reatt_Genes(x_combo)
    
      // Flatten the Reatt_Genes outputs
// Flatten Nextflow outputs first
////def flatMetas = Reatt_Genes.out.path.flatten()
////def flatDups  = Reatt_Genes.out.dup.flatten()

// Build metas channel: emits [ baseKey, folder_path ]
////def metas = flatMetas.map { p ->
////    def fullKey = p.toString().tokenize('/').find { it.startsWith('metafiles') }
////    def baseKey = fullKey.replaceAll(/(_\d+)+$/, '')  // remove numeric suffix
////    tuple(baseKey, p)                                // ✅ tuple, not nested list
////}

// Build dups channel: emits [ baseKey, dup_path ]
////def dups = flatDups.map { d ->
////    def fullKey = d.toString().tokenize('/').find { it.startsWith('dup') }
////    def baseKey = fullKey?.replace('dup', 'metafiles')
////    tuple(baseKey, d)                                // ✅ tuple, not nested list
////}
//dups.collect().ciew()
// Combine metas with their matching dup(s)
////def met_ = metas
////    .combine(dups)                       // produce all pairs
////    .filter { m -> m[0] == m[2] }     // keep only matches on baseKey
////    .map { m ->                       // 
////        def key         = m[0]
////        def folder_path = m[1]
////        def dup_path    = m[3]
////
////        // Assign CADD score
////        def cadd_score = (key == 'metafilesALL') ? 'ALL' :
////                         (key == 'metafiles20') ? '20' :
////                         (key == 'metafiles15') ? '15' : 'ALL'
////
////        tuple(folder_path, params.chromosomes, cadd_score, params.genepy_py, params.kary, dup_path)
////    }
////    .view()

// Dups as their own channel if needed
////def dup_ = dups.map { key, dup_path ->
////    def cadd_score = (key == 'metafilesALL') ? 'ALL' :
////                     (key == 'metafiles20') ? '20' :
////                     (key == 'metafiles15') ? '15' : 'ALL'
////
////    tuple(dup_path, params.chromosomes, cadd_score, params.genepy_py, params.kary, dup_path)
////}
////.view()



     //Genepy_score(dup_)
////     Genepy_score(met_)
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






