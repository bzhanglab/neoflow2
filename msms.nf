process extract_exp_info {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'

  input:
    path('manifest.txt')
  
  output:
    path "exp_*.tsv", emit: exp_ch

  """
  #!/usr/bin/env Rscript
  library(tidyverse)

  mani <- read_tsv("manifest.txt", col_types = cols(experiment = col_character()))
  exp <- mani %>% group_by(experiment)   
  experiments <- group_split(exp)
  for (i in 1:length(experiments)){
    write_tsv(experiments[[i]], paste0("exp_", i, ".tsv"), col_names=FALSE)
  }
  """
}

process extract_mgf_path {
  label 'r5_2xlarge'
  container "${params.container.ubuntu}"
  cpus 8
  memory '60 GB'

  input:
    path('exp_info.txt')
  
  output:
    tuple env(experiment_name),
          env(mgf_path),
          emit: mgf_ch

  """
  while IFS=\$'\\t' read -r -a row
  do
    experiment_name="\${row[1]}"
    mgf_path="\${row[5]}"
    break
  done < exp_info.txt
  """
}


process msms_searching{
  label 'r5_2xlarge'
  container "${params.container.neoflow}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/msms_searching/",
              mode: 'copy',
              pattern: 'exp_*/*',
              overwrite: true

  input:
    tuple val(experiment_name),  
          path(ms_file),
          path(search_db)
    path(msms_para_file)

  output:
    tuple val(experiment_name),
          path("exp_*/${res_file}"),
          emit: psm_raw_files

  script:
    if ("${params.search_engine}" == "msgf") {
        res_file = "${ms_file.baseName}.mzid"
        """
        #!/bin/sh
        java -Xmx${params.search_engine_mem}g -jar /opt/MSGFPlus.jar \
            -thread ${task.cpus} \
            -s ${ms_file} \
            -d ${search_db} \
            -conf ${msms_para_file} \
            -tda 0 \
            -o ${ms_file.baseName}.mzid
        mkdir exp_${experiment_name}
        mv *.mzid exp_${experiment_name}
        """
    }else if("${params.search_engine}" == "comet") {
        res_file = "${ms_file.baseName}_rawResults.txt"
        """
        #!/bin/sh
        /opt/comet.2018014.linux.exe -P${msms_para_file} -N${ms_file.baseName}_rawResults -D${search_db} ${ms_file}
        sed -i '1d' ${ms_file.baseName}_rawResults.txt
        sed -i '1 s/\$/\tna/' ${ms_file.baseName}_rawResults.txt
        mkdir exp_${experiment_name}
        mv *_rawResults.txt exp_${experiment_name}
        """
    }else if("${params.search_engine}" == "xtandem"){
        ms_file_name = "${ms_file.baseName}"
        xml_input = "${ms_file.baseName}_input.xml"
        res_file = "${ms_file_name}.mzid"
        """
        python3 /opt/neoflow/bin/generate_xtandem_para_xml.py ${msms_para_file} ${ms_file} ${search_db} ${ms_file_name} ${xml_input}
        ## users must provide the main search parameter file for X!Tandem search
        /opt/tandem-linux-17-02-01-4/bin/tandem.exe ${xml_input}
        ## convert xml to mzid
        java -Xmx${params.search_engine_mem}g -jar /opt/mzidlib-1.7/mzidlib-1.7.jar Tandem2mzid \
            ${ms_file_name}.xml \
            ${ms_file_name}.mzid \
            -outputFragmentation false \
            -decoyRegex XXX_ \
            -databaseFileFormatID MS:1001348 \
            -massSpecFileFormatID MS:1001062 \
            -idsStartAtZero false \
            -compress false \
            -proteinCodeRegex "\\S+"
        mkdir exp_${experiment_name}
        mv *.mzid exp_${experiment_name}
        """
  }
}


process calculate_fdr{
  label 'r5_2xlarge'
  container "${params.container.fdr_calc}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/fdr_estimation/",
              mode: 'copy',
              pattern: 'exp_*/*',
              overwrite: true

  input:
    tuple val(experiment_name),
          path(psm_raw_file),
          path(search_db)

  output:
    tuple val(experiment_name),
          path('exp_*/*_level'),
          emit: pga_result_folder
    tuple val(experiment_name),
          path('exp_*/*-rawPSMs.txt'),
          emit: raw_psm_file

  script:
    """
    mkdir -p peptide_level/global_fdr
    mkdir -p psm_level/global_fdr
    Rscript /usr/src/fdr_calc.R ./ ${search_db} ${params.pga_prefix} ${params.search_engine} ./
    mkdir exp_${experiment_name}
    mv peptide_level exp_${experiment_name}
    mv psm_level exp_${experiment_name}
    mv *-rawPSMs.txt exp_${experiment_name}
    """
}


process prepare_pepquery_input{
  label 'r5_2xlarge'
  container "${params.container.pga}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(experiment_name),
          path(pga_result_folder)

  output:
    tuple val(experiment_name),
          path("novel_peptides_psm.tsv"),
          emit: novel_psm_tsv
    tuple val(experiment_name),
          path("novel_peptides.tsv"),
          emit: novel_peptide_tsv

  """
  #!/usr/bin/env /usr/local/bin/Rscript
  library(dplyr)
  library(readr)
  library(stringr)
  is_novel_peptides=function(protein_ids){
      ids <- str_split(protein_ids,pattern=";") %>% unlist()
      is_novel <- all(str_detect(ids,pattern="^VAR"))
      return(is_novel)
  }
  
  psms <- read_tsv("peptide_level/global_fdr/pga-peptideSummary.txt")
  novel_psms <- psms %>% filter(Qvalue<=0.01) %>% 
      mutate(is_novel=sapply(protein,is_novel_peptides)) %>%
      filter(is_novel==TRUE) %>%
      filter(isdecoy==FALSE)
  ## output
  out_novel_psm_file <- "novel_peptides_psm.tsv"
  novel_psms %>% write_tsv(out_novel_psm_file)
  novel_psms %>% select(peptide) %>% distinct() %>% write_tsv("novel_peptides.tsv",col_names=FALSE)
  """
}


process run_pepquery{
  label 'r5_2xlarge'
  container "${params.container.neoflow}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/pepquery/",
              mode: 'copy',
              pattern: '*/pepquery',
              overwrite: true

  input:
    tuple val(experiment_name),
          path(ms_data),
          path(pv_refdb), 
          path(novel_peptide_tsv)

  output:
    tuple val(experiment_name),
          path("*/pepquery"),
          emit: pepquery_out

  """
  java -Xmx48g -jar /opt/pepquery-1.6.0/pepquery-1.6.0.jar \
        -ms ${ms_data} \
        -pep ${novel_peptide_tsv} \
        -db ${pv_refdb} \
        -fixMod ${params.pv_fixmod} \
        -varMod ${params.pv_varmod} \
        -cpu ${task.cpus} \
        -minScore 12 \
        -tol ${params.pv_tol} \
        -tolu ${params.pv_tolu} \
        -itol ${params.pv_itol} \
        -n 10000 \
        -um \
        -m 1 \
        -prefix neoflow \
        -o pepquery
  outdir="exp_${experiment_name}"
  mkdir \${outdir}
  mv pepquery \${outdir}
  """
}


process add_pepquery_validation {
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'
  container "${params.container.pga}"
  publishDir "${params.outdir_run}/novel_peptide_identification/",
              mode: 'copy',
              pattern: 'exp_*/*',
              overwrite: true

  input:
    tuple val(experiment_name),
          path(novel_psm_tsv),
          path(pepquery_res_folder)

  output:
    tuple val(experiment_name),
          path('*/novel_peptides_psm_pepquery.tsv'), 
          emit: novel_peptides_psm_pepquery

  """
  #!/usr/bin/env /usr/local/bin/Rscript
  library(dplyr)
  library(readr)
  library(stringr)
  psms <- read_tsv("${novel_psm_tsv}")
  psm_rank_file = "${pepquery_res_folder}/psm_rank.txt"
  if(file.exists(psm_rank_file)){
      psm_rank <- read_tsv(psm_rank_file)
      if("n_ptm" %in% names(psm_rank)){
          psm_rank <- psm_rank %>% filter(pvalue<=0.01,n_ptm==0,rank==1)
          psms\$pepquery <- ifelse(psms\$peptide %in% psm_rank\$peptide,1,0)
      }else{
         psms\$pepquery <- 0
      }
  }else{
      psms\$pepquery <- 0
  }
  dir.create(paste0("exp_", "${experiment_name}"))
  psms %>% write_tsv(file.path(paste0("exp_", "${experiment_name}"), 
                      "novel_peptides_psm_pepquery.tsv"))
  """
}



workflow msms_search {
  take:
     search_db_ch
     ref_ch

  main:
    extract_exp_info(params.manifest)
    extract_mgf_path(extract_exp_info.out.exp_ch.flatten())
    ch_1 = extract_mgf_path.out.mgf_ch
    ch_2 = search_db_ch
    search_ch = ch_1.combine(ch_2, by:0)
    msms_searching(search_ch, params.search_para_file)
    
    ch_1 = msms_searching.out.psm_raw_files
    fdr_ch = ch_1.combine(search_db_ch, by:0)
    calculate_fdr(fdr_ch)
    prepare_pepquery_input(calculate_fdr.out.pga_result_folder)

    mgf_ch = extract_mgf_path.out.mgf_ch
    pepquery_in_1 = mgf_ch.combine(ref_ch, by:0)
    novel_peptide_ch = prepare_pepquery_input.out.novel_peptide_tsv
    pepquery_in = pepquery_in_1.combine(novel_peptide_ch, by:0)
    run_pepquery(pepquery_in)

    novel_psm_ch = prepare_pepquery_input.out.novel_psm_tsv
    pepquery_res_ch = run_pepquery.out.pepquery_out
    pepquery_valid_ch = novel_psm_ch.combine(pepquery_res_ch, by:0)
    add_pepquery_validation(pepquery_valid_ch)
  
  emit:
    var_pep_file = add_pepquery_validation.out.novel_peptides_psm_pepquery
}
