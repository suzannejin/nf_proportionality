#!/usr/bin/env nextflow


/*
 * defaults parameter definitions
 */

// id
params.id = "SRP012682"

// input gtex data {id,  or already dowloaded file }
params.data = "/users/cn/sjin/projects/proportionality/data/${params.id}/rse_gene.Rdata"

// input tissue directory
params.tissue = "/users/cn/sjin/projects/proportionality/data/${params.id}/tissues/tissue10"

// if test
// params.test_size = 20000
// params.test = "--test --test_size ${params.test_size}"
params.test = [2000,5000,10000,15000,20000]

// output directory
params.outdir = "${baseDir}/results/${params.id}"



log.info """\
         PROPORTIONALITY - COMPUTE ANALYSIS FOR TISSUE SUBMATRICES"
         ======================================="
         Input dataset ID                               : ${params.id}
         Input dataset file (DATASET)                   : ${params.data}
         Input tissue file (FILE)                        : ${params.tissue}
         If test                                        : ${params.test}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()



/*
 * Managing channels
 */

// Channels containing tissue files
Channel
  .fromPath(params.tissue)
  .map { item -> [item.baseName, item] }   // Eg. [tissue1, path/to/tissue1]
  .set { ch_tissue }


/*
 * Compute proportionality for individual tissue
 */

process pr_tissue {

    // label 'process_memory'
    tag "${tissue_id}"
    publishDir "${params.outdir}/test_${test}/tissues", mode: 'copy', overwrite: true

    input:
      set val(tissue_id), \
          file(tissue_file) \
          from(ch_tissue)  
      path(data) from params.data
      each test from params.test
       
    output:
      set val(tissue_id), \
          file("${tissue_id}*.csv")

    script:
      """
      
      # print cluster working node's name
      hostname
      
      # load modules
      #module load R/3.4.1-R
      #module load OpenBLAS/0.2.15-GCC-4.9.3-2.25-LAPACK-3.6.0
      
      # print R version
      R --version
      Rscript --version
      
      # run proportionality analysis
      /usr/bin/Rscript ${baseDir}/bin/compute_tissue.R -i ${data} \
                                                       -t ${tissue_file} \
                                                       -p ${tissue_id} \
                                                       --cutoff_interval 0.01 \
                                                       --test --test_size ${test}
      
      """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
