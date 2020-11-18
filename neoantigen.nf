def extractSampleName(file_path, pattern) {
    String[] str = file_path.split('/');
    file_name = str[-1]
    result = (file_name =~ /(.*)${pattern}/)
    return result[0][1]
}


process prepare_netmhc {
  label 'r5_2xlarge'
  container "${params.container.ubuntu}"
  cpus 8
  memory '60 GB'

  input:
    path('netmhc.tgz')

  output:
    path 'netmhc', emit: netmhc_ch
  
  """
  mkdir netmhc
  tar xvfz netmhc.tgz -C netmhc --strip-components 1
  """
}


process split_file {
  label 'r5_2xlarge'
  container "${params.container.pga}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_name),
          path(var_info_file)

  output:
    tuple val(sample_name),
          path("var_info_*"),
          emit: var_info_files 

"""
  #!/usr/bin/env /usr/local/bin/Rscript
  library(dplyr)
  library(readr)
  library(parallel)
  ncpu <- detectCores()
  use_ncpu <- 0
  user_ncpu <- as.numeric("${task.cpus}")
  if(user_ncpu <= ncpu){
    use_ncpu <- user_ncpu
  }
  if(use_ncpu <= 0){
    use_ncpu <- ncpu
  }
  a <- read_tsv("${var_info_file}")
  if(use_ncpu > nrow(a)){
      ## file is small
      use_ncpu <- 1
  }
  # nlines_per_file <- ceiling(nrow(a)/use_ncpu)
  nlines_per_file <- nrow(a) %/% use_ncpu
  last_i <- 0
  for(i in 1:use_ncpu){
      i1 <- last_i + 1
      if(i < use_ncpu){    
          i2 <- i1 + nlines_per_file - 1
      }else{
          i2 <- nrow(a)
      }
      write_tsv(a[i1:i2,], paste("var_info_",i,sep=""))
      last_i <- i2
  }
  """
}


process mhc_peptide_binding_prediction {
  label 'r5_2xlarge'
  container "${params.container.binding_prediction}"
  // publishDir "${params.outdir_run}/binding_prediction/",
  //            mode: 'copy',
  //            pattern: '*/*',
  //            overwrite: true
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_name),
         path(var_db),
         path('var_info_file_??'),
         path(hla_type_file)
    path(netmhcpan_dir)

  output:
    tuple val(sample_name), 
          path("var_info_file_*_binding_prediction_result.csv"),
          emit: binding_res

  script:
  """
  chmod 755 netmhc/netMHCpan*
  chmod 755 netmhc/Linux_x86_64/bin/*
  find * -name 'var_info_file_*' | 
  parallel -j ${task.cpus} python /usr/local/bin/binding_prediction.py \
    -p {} \
    -hla_type ${hla_type_file} \
    -var_db ${var_db} \
    -var_info {} \
    -o ./ \
    -netmhcpan "${netmhcpan_dir}/netMHCpan" \
  """
}


process combine_prediction_results {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_name),
          path("*")

  output:
    tuple val(sample_name),
          file("${params.neoantigen_output_prefix}_binding_prediction_result.csv"),
          emit: res_ch

  script:
  """
  #!/usr/bin/env Rscript
  library(dplyr)
  library(readr)
  fs <- list.files(path="./",pattern="*_binding_prediction_result.csv")
  a <- lapply(fs,read.csv,stringsAsFactors=FALSE, colClasses=c("Ref"="character", "Alt"="character","AA_before"="character","AA_after"="character")) %>% bind_rows()
  ofile <- paste("${params.neoantigen_output_prefix}","_binding_prediction_result.csv",sep="")
  write_csv(a, ofile)
  """
}


/*
 * map neoepitopes to reference protiens and remove
 * neoepitopes who can map to a reference protein.
 * data preparation
 */
process prepare_data_for_mapping {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_name),
          path(mhc_binding_prediction_file)

  output:
    tuple val(sample_name),
          path("all_neoepitope.txt"),
          emit: all_neoepitope_file

  script:
  """
  #!/usr/bin/env Rscript
  library(dplyr)
  library(readr)
  a <- read_csv("${mhc_binding_prediction_file}")
  pep <- a %>% select(Neoepitope) %>% distinct()
  write_tsv(pep,"all_neoepitope.txt",col_names=FALSE)
  """
}


/*
 * map neoepitopes to reference protiens and remove
 * neoepitopes who can map to a reference protein.
 * mapping
 */
process peptide_mapping {
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'
  container "${params.container.neoflow}"

  input:
      tuple val(sample_name),
            path(all_neoepitope_file),
            path(ref_db)

  output:
      tuple val(sample_name),
            path("pep2pro.tsv"),
            emit: pep2pro

  """
  java -jar /opt/pepmap.jar -i ${all_neoepitope_file} -d ${ref_db} -o pep2pro.tsv
  """
}


/*
 * map neoepitopes to reference protiens and remove
 * neoepitopes who can map to a reference protein.
 * filtering
 */
process filtering_by_reference {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_name),
          path(pep2pro_file),
          path(mhc_binding_prediction_file)

  output:
    tuple val(sample_name), 
          path("${params.neoantigen_output_prefix}_neoepitope_filtered_by_reference.csv"),
          emit: mhc_binding_prediction_filtered_file

  """
  #!/usr/bin/env Rscript
  library(dplyr)
  library(readr)
  a <- read_csv("${mhc_binding_prediction_file}")
  pep2pro <- read_tsv("${pep2pro_file}")
  a_filter <- a %>% filter(!(Neoepitope %in% pep2pro\$peptide))
  write_csv(a_filter,"${params.neoantigen_output_prefix}_neoepitope_filtered_by_reference.csv")
  """
}


process add_variant_pep_evidence {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'
  publishDir "${params.outdir_run}/neoantigen_prediction/", mode: "copy", overwrite: true

  input:
    tuple val(sample_name),
          path(mhc_binding_prediction_filtered_file),
          path(var_pep_file),
          path(var_pep_info)

  output:
    tuple val(sample_name),
          path("*/${params.neoantigen_output_prefix}_neoepitope_filtered_by_reference_add_variant_protein_evidence.tsv"),
          emit: final_res

  """
  #!/usr/bin/env Rscript
  library(dplyr)
  library(readr)
  library(tidyr)
  a <- read.csv("${mhc_binding_prediction_filtered_file}",stringsAsFactors=FALSE) %>%
    mutate(Chr = as.character(Chr), 
            Start = as.character(Start), 
            End = as.character(End), 
            Ref = as.character(Ref),
            Alt = as.character(Alt))
  var_pep_psms <- read.delim("${var_pep_file}",stringsAsFactors=FALSE)
  var_pep_info <- read.delim("${var_pep_info}",stringsAsFactors=FALSE) %>%
    mutate(Chr = as.character(Chr))
  var_pep_pro <- var_pep_psms %>% filter(pepquery==1) %>%
    select(peptide,protein) %>% distinct() %>%
    separate_rows(protein,sep=";")
  var_pep_pro_info <- merge(var_pep_pro,var_pep_info,by.x="protein",by.y="Variant_ID") %>%
    select(peptide,Chr,Start,End,Ref,Alt) %>%
    mutate(Chr = as.character(Chr), 
            Start = as.character(Start), 
            End = as.character(End), 
            Ref = as.character(Ref),
            Alt = as.character(Alt))
  a_var <- left_join(a,var_pep_pro_info,by=c("Chr","Start","End","Ref", "Alt")) %>%
    mutate(protein_var_evidence_pep=ifelse(is.na(peptide),"-",peptide)) %>%
    mutate(peptide=NULL)
  dir.create("${sample_name}")
  a_var %>% write_tsv(file.path("${sample_name}", 
      "${params.neoantigen_output_prefix}_neoepitope_filtered_by_reference_add_variant_protein_evidence.tsv"))
  """
}



workflow neo_antigen {
  take:
     hla_type
     var_info
     var_db
     ref_ch
     var_pep_info
     var_pep_file

  main:
    prepare_netmhc(params.netmhc_file)
    // for now the type is "somatic"
    var_info.transpose()
            .map{[extractSampleName("${it[1]}", '-new-varInfo.txt'), it[1]]}
            .set{var_info_new}
    split_file(var_info_new)

    sample_exp_mapping = []
    File csvFile = new File(params.manifest)
    csvFile.eachLine { line, number ->
      if (number == 1) return
      def parts = line.split("\t")
      sample_exp_mapping.add([parts[1], parts[0]])
    }
    var_db_ch = Channel.from(sample_exp_mapping)
                    .combine(var_db)
                    .filter{it[0]==it[2]}
                    .map{[it[1],it[3]]}
    new_ch_2 = var_db_ch.combine(split_file.out.var_info_files, by:0)
    new_ch_2 = new_ch_2.combine(hla_type, by:0)
    mhc_peptide_binding_prediction(
      new_ch_2, prepare_netmhc.out.netmhc_ch
    )   
    combine_prediction_results(mhc_peptide_binding_prediction.out.binding_res)
    prepare_data_for_mapping(combine_prediction_results.out.res_ch)
    new_ref_ch = Channel.from(sample_exp_mapping)
                    .combine(ref_ch)
                    .filter{it[0]==it[2]}
                    .map{[it[1],it[3]]}
    input_ch = prepare_data_for_mapping.out.all_neoepitope_file
                 .combine(new_ref_ch, by:0)
    peptide_mapping(input_ch)
    filter_in_ch = peptide_mapping.out.pep2pro
                      .combine(combine_prediction_results.out.res_ch, by:0)
    filtering_by_reference(filter_in_ch)
    // filtering_by_reference.out.mhc_binding_prediction_filtered_file.view()
    new_var_pep_file = Channel.from(sample_exp_mapping)
                    .combine(var_pep_file)
                    .filter{it[0]==it[2]}
                    .map{[it[1],it[3]]}
    //  new_var_pep_file.view() 
     new_var_pep_info = Channel.from(sample_exp_mapping)
                    .combine(var_pep_info)
                    .filter{it[0]==it[2]}
                    .map{[it[1],it[3]]}
    //  new_var_pep_info.view() 
     new_in_ch = filtering_by_reference.out.mhc_binding_prediction_filtered_file
                    .combine(new_var_pep_file, by:0)
                    .combine(new_var_pep_info, by:0)
     add_variant_pep_evidence(new_in_ch)
}
