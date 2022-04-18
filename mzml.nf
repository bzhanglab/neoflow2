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
              pattern: '*.tar',
              overwrite: true

  input:
    tuple val(experiment), path('exp_info.txt')
  
  output:
    path "*.tar", emit: tar_ch
    path "OK", emit: ok_ch

  """
  filename="exp_info.txt"
  output_dir="exp_${experiment}"
  mkdir -p \${output_dir}
  head -n 1 \$filename | while read -r sample experiment wxs_file_name wxs_file_location mzml_files mzml_links fusion 
  do 
    IFS=';' read -r -a allUrls <<< "\$mzml_links"
    IFS=';' read -r -a allNames <<< "\$mzml_files"
    for index in "\${!allUrls[@]}"
    do
      eachUrl="\${allUrls[\$index]}"
      eachName="\${allNames[\$index]}"
      echo "Downloading \$eachUrl........."
      # wget -c "\$eachUrl" -O \${eachName} # -c allows resuming the failed or stopped downloads
      curl --connect-timeout 5  --max-time 10 --retry 5  --retry-delay 0 --retry-max-time 40 \
        --output \${eachName} \$eachUrl
      sleep 5
    done   # -c also eliminates redownload if a file already exists (with same size) in the current directory.
  done 
  mv *.gz \${output_dir}
  tar cvf ${experiment}.tar \${output_dir}
  echo "OK" > OK
  """
}


process update_manifest_file {
  label 'r5_2xlarge'
  container "${params.container.python}"
  cpus 8
  memory '60 GB'

  input:
    path('manifest.txt')
    path('OK')

  
  output:
    path "manifest_new.txt", emit: manifest_ch

  """
  #!/usr/bin/env python
  import pandas as pd

  # must start with s3, end with '/'
  s3_prefix = '${params.mzml_s3_prefix}'

  df = pd.read_csv('manifest.txt', sep='\t')
  for index, row in df.iterrows():
      exp_name = row['experiment']
      row['mzml_links'] = s3_prefix + exp_name + '.tar'

  df.to_csv('manifest_new.txt', sep='\t', index=False)
  """
}


workflow download_mzml {
  main:
    extract_exp_info(params.manifest)
    set_exp_info(extract_exp_info.out.exp_ch.flatten())
    download_mzml_files(set_exp_info.out.res_ch)
    update_manifest_file(params.manifest, download_mzml_files.out.ok_ch.collect())
  
  emit:
    manifest_new = update_manifest_file.out.manifest_ch 
}
