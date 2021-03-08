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


process download_mzml_files {
  label 'r5_2xlarge'
  container "${params.container.neoflow}"
  cpus 8
  memory '60 GB'
  publishDir "${params.mzml_s3_prefix}",
              mode: 'copy',
              pattern: '*.gz',
              overwrite: true

  input:
    tuple val(experiment), path('exp_info.txt')
  
  output:
    path("*.gz")

  """
  filename="exp_info.txt"
  # output_dir="exp_${experiment}"
  # mkdir -p \${output_dir}
  head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_uuid mzml_files mzml_links maf_file
  do 
    IFS=';' read -r -a allUrls <<< "\$mzml_links"
    IFS=';' read -r -a allNames <<< "\$mzml_files"
    for index in "\${!allUrls[@]}"
    do
      eachUrl="\${allUrls[\$index]}"
      eachName="\${allNames[\$index]}"
      echo "Downloading \$eachUrl........."
      wget -c "\$eachUrl" -O \${eachName} # -c allows resuming the failed or stopped downloads
      sleep 5
    done   # -c also eliminates redownload if a file already exists (with same size) in the current directory.
  done 
  """
}


workflow mzml {
  main:
    extract_exp_info(params.manifest)
    set_exp_info(extract_exp_info.out.exp_ch.flatten())
    download_mzml_files(set_exp_info.out.res_ch)
}
