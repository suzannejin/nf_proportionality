#!/usr/bin/env nextflow


/*
 * defaults parameter definitions
 */

// id
params.id = "SRP012682"

// input gtex data {id,  or already dowloaded file }
params.data = "/users/cn/sjin/projects/proportionality/data/${params.id}/rse_gene.Rdata"

// input donor directory
params.donor = "/users/cn/sjin/projects/proportionality/data/${params.id}/donors/donor1?"
// params.donor = "/users/cn/sjin/projects/proportionality/data/${params.id}/donors/donor1"

// if test
// params.test_size = 20000
// params.test = "--test --test_size ${params.test_size}"
params.test = ""

// only results
params.only_results = "--only_results"
// params.only_results = ""

// output directory
params.outdir = "${baseDir}/results/${params.id}"



log.info """\
         PROPORTIONALITY - COMPUTE ANALYSIS FOR DONORS SUBMATRICES"
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input donor file (FILE)                        : ${params.donor}
         If test                                        : ${params.test}
         If only results                                : ${params.only_results}
         Output directory (DIRECTORY)                   : ${params.outdir}
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
      path(data) from params.data
      val(chunk_size) from params.chunk_size
      val(test) from params.test
      val(only_results) from params.only_results
       
    output:
      set val(donor_id), \
          file("${donor_id}*.csv")

    script:
      """
      
      # print cluster working node's name
      hostname
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      Rscript ${baseDir}/bin/compute_donor.R -i ${data} \
                                             -d ${donor_file} \
                                             -p ${donor_id} \
                                             --cutoff_interval 0.01 \
                                             --chunk_size ${chunk_size}
                                             ${test} ${only_results}
      
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
