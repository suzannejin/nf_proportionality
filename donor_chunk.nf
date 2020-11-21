#!/usr/bin/env nextflow


/*
 * defaults parameter definitions
 */

// id
params.id = "SRP012682"

// input gtex data {id,  or already dowloaded file }
params.data = "/users/cn/sjin/projects/proportionality/data/${params.id}/rse_gene.Rdata"

// input donor directory
params.donor = "/users/cn/sjin/projects/proportionality/data/${params.id}/donors/donor56"

// chunk size
params.chunk_size = 5000

// if test
// params.test_size = 20000
// params.test = "--test --test_size ${params.test_size}"
params.test = ""

// output directory
params.outdir = "${baseDir}/results_chunk3/${params.id}"

// if check output file exists
params.check = true



log.info """\
         PROPORTIONALITY - COMPUTE ANALYSIS FOR DONORS SUBMATRICES"
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input donor file (FILE)                        : ${params.donor}
         Chunk size (int)                               : ${params.chunk_size}
         If test                                        : ${params.test}
         Output directory (DIRECTORY)                   : ${params.outdir}
         If check output exists (bool)                  : ${params.check}
         """
         .stripIndent()



/*
 * Managing channels
 */

// Channels containing donor files
Channel
  .fromPath(params.donor)
  .map { item -> [item.baseName, item] }   // Eg. [donor3, path/to/donor3]
  .set { ch_donor }


/*
 * Compute proportionality for individual donors
 */

process pr_donor {

    // label 'process_memory'
    tag "${donor_id}"
    publishDir "${params.outdir}/donors", mode: 'copy', overwrite: true

    input:
      set val(donor_id), \
          file(donor_file) \
          from(ch_donor)  
      val(chunk_size) from params.chunk_size
      path(data) from params.data
      val(test) from params.test
       
    output:
      set val(donor_id), \
          file("${donor_id}*.csv")

    when:
      if(params.check){
        if(!file("${params.outdir}/donors/${donor_id}_results.csv").exists()){true}else{false}
      }else{
        true
      }

    script:
      """
      
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/compute_donor_chunk.R -i ${data} \
                                                   -d ${donor_file} \
                                                   -p ${donor_id} \
                                                   --cutoff_interval 0.01 \
                                                   --ncores ${task.cpus} \
                                                   --chunk_size ${chunk_size} \
                                                   ${test} 
      
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
