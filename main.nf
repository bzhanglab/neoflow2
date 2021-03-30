#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

assert params.run_version
assert params.manifest
assert params.maf
assert params.fusion_file
// the first step is to download the mzml files to s3 
// this is where downloaded files will be stored first 
assert params.mzml_s3_prefix
// these have default values
assert params.search_engine
assert params.search_para_file

params.outdir_run = "${params.outdir}/${params.run_version}"

def helpMessage() {
    log.info neoflowHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile docker
    or 
    nextflow run main.nf -profile awsbatch -bucket-dir s3://zhanglab-nextflow-workdir/workdir
    Mandatory arguments:
      -profile                      Configuration profile to use.
                                    Available: docker, awsbatch.
    """.stripIndent()
}

header = """==============================================
                  ___ __               ___ 
  ___  ___  ___   / _// /___  _    __  |_  |
 / _ \\/ -_)/ _ \\ / _// // _ \\| |/|/ / / __/ 
/_//_/\\__/ \\___//_/ /_/ \\___/|__,__/ /____/ 

==============================================
"""

def neoflowHeader() {
    return header.stripIndent()
}

include { download_mzml } from './mzml'
include { hla_typing } from './hla_typing'
include { database_construction } from './db_construction'
include { msms_search } from './msms'
include { neo_antigen } from './neoantigen'


def check_input_files() {
   log.info 'checking input files....'
   manifest_file = file("${params.manifest}")
   if (!manifest_file.isFile()) {
     log.info 'manifest file not exist'
     return false
   }
   fusion_tgz = file("${params.fusion_file}")
   if (!fusion_tgz.isFile()) {
     log.info 'fusion tgz file not exist'
     return false
   }
   fusion_tgz.copyTo('./fusion_file.tgz')
   myDir = file('./fusion')
   result = myDir.mkdir()
   cmd = 'tar xvfz fusion_file.tgz -C ./fusion --strip-components=1'
   def proc = cmd.execute()
   proc.waitForOrKill(5000)

   count = 0
   allLines  = manifest_file.readLines()
   status = true
   for( line : allLines ) {
      count = count + 1
      if (count == 1) continue
      String[] str
      str = line.split('\t')
      fusion_path = str[6]
      cur_file = file(fusion_path)
      if (!cur_file.isFile()) {
         log.error "${fusion_path} does not exist"
         if (status) {
           status = false
         }
      }
    }

   // clean up
   cmd = '/bin/rm -rf ./fusion ./fusion_file.tgz'
   cmd.execute()
   log.info 'file checking done...'
   return status
}


workflow neoflow2_sub {
    // mzml source is from PDC,
    // first thing to do is download the data to s3,
    // then the manifest file will be updated with the new s3 URI
    download_mzml()
    hla_typing(download_mzml.out.manifest_new)
    database_construction(download_mzml.out.manifest_new)
    msms_search(
      download_mzml.out.manifest_new,
      database_construction.out.search_db_ch,
      database_construction.out.ref_ch
    )
    neo_antigen(
      hla_typing.out.hla_typing_out,
      database_construction.out.sample_varinfo_ch,
      database_construction.out.var_db_ch,
      database_construction.out.ref_ch,
      database_construction.out.var_pep_info,
      msms_search.out.var_pep_file
    )
}


// ====== main workflow ===========
workflow {
  log.info neoflowHeader()
  if (check_input_files()) {
     neoflow2_sub()
  } 
}


report_txt = """$header
 name = ${workflow.manifest.name}
 author = ${workflow.manifest.author}
 homePage = ${workflow.manifest.homePage}
 description = ${workflow.manifest.description}
 mainScript = ${workflow.manifest.mainScript}
 nextflowVersion = ${workflow.manifest.nextflowVersion}
 version = ${workflow.manifest.version}
""".stripIndent()


workflow.onComplete {
  def output_d = new File("results/pipeline_info/")
  if (!output_d.exists()) {
      output_d.mkdirs()
  }
  println "Pipeline completed at: $workflow.complete"
  def output_tf = new File(output_d, "pipeline_report.txt")
  output_tf.withWriter { w -> w << report_txt }

}