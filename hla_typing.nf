process generate_id_files {
  container "${params.container.r_tidyverse}"
  cpus 1
  memory '4 GB'
  executor 'local'

  input:
    path('manifest.txt')
    val start
    val end

   output:
    path 'case_*.tsv'

  script:
  """
   #!/usr/bin/env Rscript
   library(tidyverse)

  # old_manifest <- read_tsv("manifest.txt")
  # mainfest <- old_manifest
  # total_n_rows <- nrow(manifest)
  # # rename sample names with unique row id R-XXXXX
  # row_ids <- paste0('R-', sprintf('%0.5d', 1:total_n_rows))
  # manifest[1] <- row_ids
  # old_manifest %>% 
  #   add_column(row_id = row_ids, .before = 'sample') %>%
  #  write_tsv('new_manifest.tsv')
       

   manifest <- read_tsv("manifest.txt")
   case_start <- ifelse($start == -1, 1, $start)
   case_end <- ifelse($end == -1, nrow(manifest), $end)
   stopifnot(case_start >= 1, case_start <= nrow(manifest))
   stopifnot(case_end >= 1, case_end <= nrow(manifest))
   stopifnot(case_start <= case_end)

   # write_tsv(all_case_tbl, "all_case_info.tsv")
   for (i in seq(case_start, case_end)) {
     cur_row <- manifest %>% slice(i)
     write_tsv(cur_row, paste0('case_', i, '.tsv'), col_name = FALSE)
   }
   """
}


process get_sample_id {
  container "${params.container.ubuntu}"
  cpus 1
  memory '4 GB'
  executor 'local'

  input:
    path(id_file)
  
  output:
    tuple env(sample_id), path("location.txt"), emit: res_ch

   """
   while IFS=\$'\\t' read -r -a myid
   do
      sample_id="\${myid[0]}"
      location="\${myid[3]}" 
   done < "${id_file}"
   echo \${location} > location.txt
   # uuid=\$(echo "\${uuid}" | sed -e 's/?/\\\\?/g' -e 's/=/\\\\=/g' -e 's/&/\\\\&/g' -e 's/%/\\\\%/g')
   """
}

process get_path_string {
  container "${params.container.ubuntu}"
  cpus 1
  memory '4 GB'
  executor 'local'

  input:
    tuple val(sample_id), path(path_file)
  
  output:
    tuple val(sample_id), env(path),  emit: res_ch

   """
   path=\$(cat ${path_file})
   """
}


process download_files_path {
  container "${params.container.neoflow}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), path(file_path)

 output:
   tuple val(sample_id), path(file_path), emit: res_ch

  """
  true
  """
}


process download_files_uuid {
  container "${params.container.gdc_client}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    path(token)
    tuple val(sample_id), path(uuid_file)

 output:
   tuple val(sample_id),
         path('bam/*.bam'),
         emit: res_ch

  """
    uuid=`cat ${uuid_file}`
    gdc-client download "\${uuid}" -t "${token}" -n ${task.cpus} --retry-amount 5 --wait-time 30
    mv "\${uuid}" bam
  """
}


process download_files_url_bam {
  container "${params.container.neoflow}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), path(url_file)

 output:
   tuple val(sample_id),
         path('bam/*.bam'),
         emit: res_ch

  """
    mkdir bam
    cd bam
    url=`cat ../${url_file}`
    curl  --retry 5 \
    --retry-delay 5 --retry-max-time 0 \
    \"\${url}\" --output ${sample_id}.bam
  """  
}


// cram crai
process download_files_url_cram {
  container "${params.container.neoflow}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), path(url_file)

 output:
   tuple val(sample_id),
         path('bam/*.cram'),
         emit: res_ch

  """
    mkdir bam
    cd bam
    url=`cat ../${url_file}`
    curl --retry 5 \
    --retry-delay 5 --retry-max-time 0 \
    \"\${url}\" --output ${sample_id}.cram
  """  
}


process bam_to_fastq {
  label 'r5_xlarge'
  container  "${params.container.samtools}"
  cpus 4
  memory '30 GB'

  input:
    tuple val(sample_id),
          path('wxs.bam')


  output:
    tuple val(sample_id), 
          path('r*{1,2}.fastq.gz'),
          emit: res_ch


  """
  samtools fastq  -@ ${task.cpus} -1 r1.fastq.gz -2 r2.fastq.gz wxs.bam
  """
}


process cram_to_fastq {
  label 'r5_xlarge'
  container  "${params.container.samtools}"
  cpus 4
  memory '30 GB'

  input:
    tuple val(sample_id),
          path('wxs.cram')

  output:
    tuple val(sample_id), 
          path('r*{1,2}.fastq.gz'),
          emit: res_ch


  """
  samtools fastq  -@ ${task.cpus} -1 r1.fastq.gz -2 r2.fastq.gz wxs.cram
  """

}


process reads_mapping {
  label 'c5a_8xlarge'
  container  "${params.container.bwa}"
  cpus 32
  memory '60 GB'

  input:
    tuple val(sample_id), path('*')
    val(hla_ref_prefix)
    path('*')


  output:
     tuple val(sample_id), path('mapped_{1,2}.fastq'), emit: res_ch

  script:
    """
    bwa mem -t ${task.cpus} -M ${hla_ref_prefix} r1.fastq.gz | samtools view -@ ${task.cpus} -bS - > mapped_1.bam
    rm -rf r1.fastq.gz
    samtools view -@ ${task.cpus} -h -F 4 -b1 mapped_1.bam > mapped_11.bam
    rm -rf mapped_1.bam
    samtools fastq -@ ${task.cpus} mapped_11.bam > mapped_1.fastq
    rm -rf mapped_11.bam
    
    bwa mem -t ${task.cpus} -M ${hla_ref_prefix} r2.fastq.gz | samtools view -@ ${task.cpus} -bS - > mapped_2.bam
    rm -rf r2.fastq.gz
    samtools view -@ ${task.cpus} -h -F 4 -b1 mapped_2.bam > mapped_22.bam
    rm -rf mapped_2.bam
    samtools fastq -@ ${task.cpus} mapped_22.bam > mapped_2.fastq
    rm -rf mapped_22.bam
    """
}


process run_optitype {
  label 'r5_2xlarge'
  container  "${params.container.optitype}"
  publishDir "${params.outdir_run}/hla_type/", 
              pattern: "optitype_results/${sample_id}",
              mode: 'copy', 
              overwrite: true
  cpus 8
  memory '60 GB'

  input:
      tuple val(sample_id), path('*')

  output:
     path "optitype_results/${sample_id}"
     tuple val(sample_id),
           path("optitype_results/${sample_id}/*.tsv"),
           emit: res_ch

  script:
    """
    python /opt/OptiType/OptiTypePipeline.py \
      --input mapped_1.fastq mapped_2.fastq \
      --prefix ${sample_id} \
      --outdir optitype_results/${sample_id} \
      --${params.seqtype}
    """
}


workflow hla_typing {
    take:
      manifest_new

    main:
      generate_id_files(manifest_new, params.start, params.end)
      get_sample_id(generate_id_files.out.flatten())
      if (params.bam_source == 'uuid') {
         download_files_uuid(params.gdc_token, get_sample_id.out.res_ch)
         fastq_out = bam_to_fastq(download_files_uuid.out.res_ch)
      } else if (params.bam_source == 'url'){ 
        if (params.bam_type == 'bam') {   // url + bam
          download_files_url_bam(get_sample_id.out.res_ch)
          fastq_out = bam_to_fastq(download_files_url_bam.out.res_ch)
        }
        else {  // url + cram
          download_files_url_cram(get_sample_id.out.res_ch)
          fastq_out = cram_to_fastq(download_files_url_cram.out.res_ch)
        }
      } else {  // path
        get_path_string(get_sample_id.out.res_ch)
        download_files_path(get_path_string.out.res_ch)
        if (params.bam_type == 'bam') { // path + bam
          fastq_out = bam_to_fastq(download_files_path.out.res_ch)
        } else {  // path + cram
          fastq_out = cram_to_fastq(download_files_path.out.res_ch)
        }
      }
      reads_mapping(
        fastq_out.res_ch,
        params.hla_ref_prefix,
        Channel.fromPath(params.hla_ref).collect()
      )
     run_optitype(reads_mapping.out.res_ch)

     emit:
       hla_typing_out = run_optitype.out.res_ch
}
