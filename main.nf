#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

assert params.run_version
assert params.manifest  
assert params.bam_source  // 'uuid' or 'url'
assert params.bam_type   // 'bam' or 'cram'

if(!params.hlatyping) {
  assert params.maf
  assert params.fusion_file
}
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
      --run_version
      --manifest
      --mzml_s3_prefix

    Optional arguments:
      --maf
      --fusion_file 
      --database
      --hlaptyping
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
   if (params.hlatyping) {
     return true
   }
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
  exp_sample_mapping = [:]
  exp_names = []
  sample_names = []
  
  csvFile = file(params.manifest)
  allLines  = csvFile.readLines()
  number = 0
  for( line : allLines ) {
    number++
    if (number == 1) continue
    def parts = line.split("\t")
    if (exp_sample_mapping.containsKey(parts[1])){
       exp_sample_mapping[parts[1]].add(parts[0])
    } else{
       exp_sample_mapping[parts[1]] = [parts[0]]
    }
    exp_names.add(parts[1])
    sample_names.add(parts[0])
  }
  sample_names = sample_names.unique().toSorted()
  exp_names = exp_names.unique().toSorted()

  // mzml source is from PDC or AWS S3 
  // the provided are urls for the files
  // first thing to do is download the data to s3,
  // then the manifest file will be updated with the new s3 URI
  download_mzml()
  if (!params.hlatyping) {
    hla_typing(download_mzml.out.manifest_new)
    hla_typing_out_ch = hla_typing.out.hla_typing_out
  } else {
    // create hlatyping related channel if use provides hla typing results
    hlatyping_input = params.hlatyping.replaceAll(/\/$/,"")
    hlatyping_ch_lst = []
    sample_names.each{
      sample_name -> 
      hlatyping_ch_lst.add(["${sample_name}", 
        "${hlatyping_input}/optitype_results/${sample_name}/${sample_name}_result.tsv"])
    }
    hla_typing_out_ch = Channel.fromList(hlatyping_ch_lst)
  }
  if (!params.database) {
    database_construction(download_mzml.out.manifest_new)
    search_db_ch = database_construction.out.search_db_ch
    ref_ch = database_construction.out.ref_ch
    sample_varinfo_ch = database_construction.out.sample_varinfo_ch
    var_db_ch = database_construction.out.var_db_ch
    var_pep_info_ch = database_construction.out.var_pep_info
  } else {
    database_input = params.database.replaceAll(/\/$/,"")
    search_db_ch_lst = []
    exp_names.each{
      exp_name -> 
      search_db_ch_lst.add(["${exp_name}", 
        "${database_input}/exp_${exp_name}/${exp_name}_target_decoy.fasta"])
    }
    search_db_ch = Channel.fromList(search_db_ch_lst)

    ref_lst = []
    exp_names.each{
      exp_name -> 
      ref_lst.add(["${exp_name}", 
        "${database_input}/exp_${exp_name}/ref.fasta"])
    }
    ref_ch = Channel.fromList(ref_lst)
    
    sample_varinfo_lst = []
    exp_names.each{
      exp_name -> 
      tmp_lst = []
      exp_sample_mapping[exp_name].each{
        sample_name -> 
        tmp_lst.add("${database_input}/exp_${exp_name}/sample_varinfo/${sample_name}-new-varInfo.txt")
      }
      sample_varinfo_lst.add(["${exp_name}", tmp_lst])
    }
    sample_varinfo_ch = Channel.fromList(sample_varinfo_lst)
    var_db_ch_lst = []
    exp_names.each{
      exp_name -> 
      var_db_ch_lst.add(["${exp_name}", 
        "${database_input}/exp_${exp_name}/${exp_name}-anno_varInfor_with_fusion.fasta"])
    }
    var_db_ch = Channel.fromList(var_db_ch_lst)
    var_pep_info_lst = []
    exp_names.each{
      exp_name -> 
      var_pep_info_lst.add(["${exp_name}", 
         "${database_input}/exp_${exp_name}/experiment_varinfo/${exp_name}-anno_varInfor_with_fusion.txt"])
    }
    var_pep_info_ch = Channel.fromList(var_pep_info_lst)
  }

  msms_search(
    download_mzml.out.manifest_new,
    search_db_ch,
    ref_ch
  )
  neo_antigen(
    hla_typing_out_ch,
    sample_varinfo_ch,
    var_db_ch,
    ref_ch,
    var_pep_info_ch,
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