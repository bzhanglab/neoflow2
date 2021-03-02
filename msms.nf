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
    exp_name <- unique(experiments[[i]][["experiment"]])
    write_tsv(experiments[[i]], paste0("exp_", exp_name, ".tsv"), col_names=FALSE)
  }
  """
}

process set_exp_info {
  label 'r5_2xlarge'
  container "${params.container.ubuntu}"
  cpus 8
  memory '60 GB'

  input:
     path(exp_info_file)
  
  output:
     tuple env(exp_name), path(exp_info_file), emit: res_ch
     
  """
  exp_info=${exp_info_file}
  exp_name=\${exp_info#exp_}
  exp_name=\${exp_name%.tsv}
  """
}


process extract_mzml_path {
  label 'r5_2xlarge'
  container "${params.container.neoflow}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/msms_searching/",
              mode: 'copy',
              pattern: 'exp_*/*.mzid',
              overwrite: true

  input:
    tuple val(experiment), path('exp_info.txt'), path('ref.fasta')
    path(msms_para_file)
  
  output:
    tuple val(experiment),
          path("exp_*/*.mzML"),
          emit: mzml_ch
    tuple val(experiment),
          path("exp_*/*.mzid"),
          emit: mzid_ch

  """
  filename="exp_info.txt"
  output_dir="exp_${experiment}"
  mkdir -p \${output_dir}
  head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_uuid mzml_files mzml_links maf_file
  do 
    IFS=';' read -r -a allUrls <<< "\$mzml_links"
    IFS=';' read -r -a allNames <<< "\$mzml_files"
    for index in "\${!allUrls[@]}"
    do
      eachUrl="\${allUrls[\$index]}"
      eachName="\${allNames[\$index]}"
      echo "Downloading \$eachUrl........."
      wget -c "\$eachUrl" -O \${output_dir}/\${eachName} # -c allows resuming the failed or stopped downloads
      gunzip \${output_dir}/\${eachName}
      msmsName=\${output_dir}/\${eachName}
      msmsName=\${msmsName/.gz/}
      java -Xmx${params.search_engine_mem}g -jar /opt/MSGFPlus.jar \
        -thread ${task.cpus}\
        -s \$msmsName \
        -d ref.fasta \
        -tda 0 \
        -o \${msmsName}.mzid \
        -conf ${msms_para_file}
    done   # -c also eliminates redownload if a file already exists (with same size) in the current directory.
  done 
  """
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
              pattern: 'exp_*/pepquery',
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
  mkdir -p pepquery_index
  java -cp /opt/pepquery-1.6.2/pepquery-1.6.2.jar main.java.index.BuildMSlibrary \
       -i ./  -o pepquery_index
  mkdir -p pepquery
  java -Xmx48g -jar /opt/pepquery-1.6.2/pepquery-1.6.2.jar \
        -ms pepquery_index \
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
    set_exp_info(extract_exp_info.out.exp_ch.flatten())
    combined_ch = set_exp_info.out.res_ch.combine(search_db_ch, by:0)
    extract_mzml_path(combined_ch, params.search_para_file)
    ch_1 = extract_mzml_path.out.mzid_ch
    fdr_ch = ch_1.combine(search_db_ch, by:0)
    calculate_fdr(fdr_ch)
    prepare_pepquery_input(calculate_fdr.out.pga_result_folder)

    mzml_ch = extract_mzml_path.out.mzml_ch
    pepquery_in_1 = mzml_ch.combine(ref_ch, by:0)
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
