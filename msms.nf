process extract_exp_info {
  container "${params.container.r_tidyverse}"
  cpus 1
  memory '4 GB'
  executor 'local'

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
  container "${params.container.ubuntu}"
  cpus 1
  memory '4 GB'
  executor 'local'

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


process extract_s3_path {
  container "${params.container.python}"
  cpus 1
  memory '4 GB'
  executor 'local'

  input:
    tuple val(experiment), path('exp_info.txt')

  output:
    tuple val(experiment), env(s3_uri), emit: s3_ch

  """
  # only read the first line
   while IFS=\$'\\t' read -r -a myline
   do
      s3_uri="\${myline[5]}"
      break
   done < exp_info.txt
  """

}


process peptide_identification {
  label 'r5_2xlarge'
  container "${params.container.neoflow}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/msms_searching/",
              mode: 'copy',
              pattern: 'exp_*/*.mzid',
              overwrite: true

  input:
    tuple val(experiment), path('exp_info.txt'), path(mzml_tar), path('ref.fasta')
    path(msms_para_file)
  
  output:
    tuple val(experiment),
          path("exp_*/*.mzML"),
          emit: mzml_ch
    tuple val(experiment),
          path("exp_*/*${output_pattern}"),
          emit: mzid_ch

  script:
      if (params.search_engine == 'msgf') {
          output_pattern = ".mzid"
          """
          # mzml_tar contains all mzml file
          filename="exp_info.txt"
          output_dir="exp_${experiment}"
          mkdir -p \${output_dir}
          ARCHIVE=${mzml_tar}
          CONTENTS=\$(tar -tf "\$ARCHIVE")
          # Check if the first entry in the archive is a directory
          FIRST_ENTRY=\$(echo "\$CONTENTS" | head -n 1)
          if [[ "\$FIRST_ENTRY" == */ ]]; then
            tar -xvf "\$ARCHIVE" --strip-components 1
          else
            tar -xvf "\$ARCHIVE"
          fi
          mv *.gz \${output_dir}
          head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_location mzml_files mzml_path fusion
          do 
            IFS=';' read -r -a allNames <<< "\$mzml_files"
            for index in "\${!allNames[@]}"
            do
              eachName="\${allNames[\$index]}"
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
            done   
          done 
          """
      } else if (params.search_engine == 'comet'){
        output_pattern = "_rawResults.txt"
        """
          # mzml_tar contains all mzml file
          filename="exp_info.txt"
          output_dir="exp_${experiment}"
          mkdir -p \${output_dir}
          ARCHIVE=${mzml_tar}
          CONTENTS=\$(tar -tf "\$ARCHIVE")
          # Check if the first entry in the archive is a directory
          FIRST_ENTRY=\$(echo "\$CONTENTS" | head -n 1)
          if [[ "\$FIRST_ENTRY" == */ ]]; then
            tar -xvf "\$ARCHIVE" --strip-components 1
          else
            tar -xvf "\$ARCHIVE"
          fi
          mv *.gz \${output_dir}
          head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_location mzml_files mzml_path fusion
          do 
            IFS=';' read -r -a allNames <<< "\$mzml_files"
            for index in "\${!allNames[@]}"
            do
              eachName="\${allNames[\$index]}"
              gunzip \${output_dir}/\${eachName}
              msmsName=\${output_dir}/\${eachName}
              msmsName=\${msmsName/.gz/}
              /opt/comet.2018014.linux.exe -P ${msms_para_file} -N \$msmsName_rawResults -D ref.fasta \${msmsName}
              sed -i '1d' \${msmsName}_rawResults.txt
              sed -i '1 s/\$/\tna/' \${msmsName}_rawResults.txt
            done
          done
        """
      } else if (params.search_engine == 'xtandem'){
        output_pattern = ".mzid"
        """
        # mzml_tar contains all mzml file
        filename="exp_info.txt"
        output_dir="exp_${experiment}"
        mkdir -p \${output_dir}
        ARCHIVE=${mzml_tar}
        CONTENTS=\$(tar -tf "\$ARCHIVE")
        # Check if the first entry in the archive is a directory
        FIRST_ENTRY=\$(echo "\$CONTENTS" | head -n 1)
        if [[ "\$FIRST_ENTRY" == */ ]]; then
          tar -xvf "\$ARCHIVE" --strip-components 1
        else
          tar -xvf "\$ARCHIVE"
        fi
        mv *.gz \${output_dir}
        head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_location mzml_files mzml_path fusion
        do 
          IFS=';' read -r -a allNames <<< "\$mzml_files"
          for index in "\${!allNames[@]}"
          do
            eachName="\${allNames[\$index]}"
            gunzip \${output_dir}/\${eachName}
            msmsName=\${output_dir}/\${eachName}
            msmsName=\${msmsName/.gz/}
            xml_input=\${msmsName}_input.xml
            python3 /opt/neoflow/bin/generate_xtandem_para_xml.py ${msms_para_file} \${msmsName} ref.fasta \${msmsName} \${xml_input}
            ## users must provide the main search parameter file for X!Tandem search
            /opt/tandem-linux-17-02-01-4/bin/tandem.exe \${xml_input}
            ## convert xml to mzid
            java -Xmx${params.search_engine_mem}g -jar /opt/mzidlib-1.7/mzidlib-1.7.jar Tandem2mzid \
              \${msmsName}.xml \
              \${msmsName}.mzid \
              -outputFragmentation false \
              -decoyRegex XXX_ \
              -databaseFileFormatID MS:1001348 \
              -massSpecFileFormatID MS:1001062 \
              -idsStartAtZero false \
              -compress false \
              -proteinCodeRegex "\\S+"
            done
          done
        """
      } else {
        """
        # exit with 1 if the search engine is not supported
        exit 1
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
  container "${params.container.pga}"
  cpus 1
  memory '4 GB'
  executor 'local'

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
  cpus 1
  memory '4 GB'
  executor 'local'
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
     manifest_new
     search_db_ch
     ref_ch

  main:
    extract_exp_info(manifest_new)
    set_exp_info(extract_exp_info.out.exp_ch.flatten())
    extract_s3_path(set_exp_info.out.res_ch)
    combined_ch_1 = set_exp_info.out.res_ch.combine(extract_s3_path.out.s3_ch, by:0)
    combined_ch = combined_ch_1.combine(search_db_ch, by:0)
    peptide_identification(combined_ch, params.search_para_file)
    ch_1 = peptide_identification.out.mzid_ch
    fdr_ch = ch_1.combine(search_db_ch, by:0)
    calculate_fdr(fdr_ch)
    prepare_pepquery_input(calculate_fdr.out.pga_result_folder)

    mzml_ch = peptide_identification.out.mzml_ch
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
