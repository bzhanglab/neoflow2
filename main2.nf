#!/usr/bin/env nextflow

// database construction and hla_typing is done 
// the results are provided 
// only perform msms_search and neo_antigen subworkflows
nextflow.enable.dsl = 2

assert params.run_version
assert params.search_engine
assert params.search_para_file
assert params.database
assert params.hlatyping

params.outdir_run = "${params.outdir}/${params.run_version}"
database_input = params.database.replaceAll(/\/$/,"")
hlatyping_input = params.hlatyping.replaceAll(/\/$/,"")

def helpMessage() {
    log.info neoflowHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile docker --database /path/to/database/results --hlatyping /path/to/hlatyping/results
    or 
    nextflow run main.nf -profile awsbatch -bucket-dir s3://zhanglab-nextflow-workdir/workdir --database /path/to/database/results --hlatyping /path/to/hlatyping/results
    Mandatory arguments:
      -profile                      Configuration profile to use.
                                    Available: docker, awsbatch.
      --database
      --hlaptyping
      --run_version
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

include { msms_search } from './msms'
include { neo_antigen } from './neoantigen'

// ====== main workflow ===========
workflow {
  log.info neoflowHeader()
  // create database construction results channel
  // create hla typing results channel
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

  



  File csvFile = new File(params.manifest)
  csvFile.eachLine { line, number ->
    if (number == 1) return
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

  hlatyping_ch_lst = []
  sample_names.each{
    sample_name -> 
    hlatyping_ch_lst.add(["${sample_name}", 
      "${hlatyping_input}/optitype_results/${sample_name}/${sample_name}_result.tsv"])
  }
  hla_typing_out = Channel.fromList(hlatyping_ch_lst)

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
      "${database_input}/exp_${exp_name}/${exp_name}_anno-var.fasta"])
  }
  var_db_ch = Channel.fromList(var_db_ch_lst)

  var_pep_info_lst = []
  exp_names.each{
    exp_name -> 
    var_pep_info_lst.add(["${exp_name}", 
      "${database_input}/exp_${exp_name}/experiment_varinfo/${exp_name}_anno-varInfo.txt"])
  }
  var_pep_info = Channel.fromList(var_pep_info_lst)

  msms_search(
    search_db_ch,
    ref_ch
  )
   neo_antigen(
     hla_typing_out,
     sample_varinfo_ch,
     var_db_ch,
     ref_ch,
     var_pep_info,
     msms_search.out.var_pep_file
   )
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