process generate_id_files {
  label 'r5_2xlarge'
  container "${params.container.r_tidyverse}"
  cpus 8
  memory '60 GB'

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
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    path(id_file)
  
  output:
    tuple env(sample_id), env(uuid), emit: res_ch

   """
   while IFS=\$'\\t' read -r -a myid
   do
      sample_id="\${myid[0]}"
      uuid="\${myid[3]}" 
   done < "${id_file}"
   """
}


process download_files_uuid {
  container "${params.container.gdc_client}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    path(token)
    tuple val(sample_id), val(location)

 output:
   tuple val(sample_id),
         path('bam/*.bam'),
         emit: res_ch

  """
    gdc-client download "${location}" -t "${token}" -n ${task.cpus}
    mv "${location}" bam
  """
}


process download_files_url_bam {
  container "${params.container.ubuntu}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), val(location)

 output:
   tuple val(sample_id),
         path('bam/*.bam'),
         emit: res_ch

  """
    mkdir bam
    cd bam
    curl --connect-timeout 5 \
    --max-time 10 --retry 5 \
    --retry-delay 0 --retry-max-time 40 \
    \${location} --output \${sample_id}.bam
  """  
}


// cram crai
process download_files_url_cram {
  container "${params.container.ubuntu}"
  label 'r5_2xlarge'
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), val(location)

 output:
   tuple val(sample_id),
         path('bam/*.cram'),
         emit: res_ch

  """
    mkdir bam
    cd bam
    curl --connect-timeout 5 \
    --max-time 10 --retry 5 \
    --retry-delay 0 --retry-max-time 40 \
    \${location} --output \${sample_id}.cram
  """  
}



process bam_to_fastq {
  label 'r5_2xlarge_500g'
  container  "${params.container.samtools}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id),
          path('wxs.bam')


  output:
    tuple val(sample_id), 
          path('r*{1,2}.fastq'),
          emit: res_ch


  """
  samtools fastq  -@ ${task.cpus} -1 r1.fastq -2 r2.fastq wxs.bam
  """
}


process cram_to_fastq {
  label 'r5_2xlarge_500g'
  container  "${params.container.samtools}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id),
          path('wxs.cram')

  output:
    tuple val(sample_id), 
          path('r*{1,2}.fastq'),
          emit: res_ch


  """
  samtools fastq  -@ ${task.cpus} -1 r1.fastq -2 r2.fastq wxs.cram
  """

}


process reads_mapping {
  label 'r5_2xlarge_500g'
  container  "${params.container.bwa}"
  cpus 8
  memory '60 GB'

  input:
    tuple val(sample_id), path('*')
    val(hla_ref_prefix)
    path('*')


  output:
     tuple val(sample_id), 
          path('mapped_{1,2}.sam'),
          emit: res_ch

  script:
    """
    bwa mem -t ${task.cpus} -M ${hla_ref_prefix} r1.fastq > mapped_1.sam
    rm -rf r1.fastq
    bwa mem -t ${task.cpus} -M ${hla_ref_prefix} r2.fastq > mapped_2.sam
    rm -rf r2.fastq
    """
}


process run_samtools{
  label 'r5_2xlarge_500g'
  container  "${params.container.samtools}"
  cpus 8
  memory '60 GB'

  input:
    path ('*')
    tuple val(sample_id), path('*') 

  output:
    tuple val(sample_id), path('mapped_{1,2}.fastq'), emit: res_ch

  script:
    """
    samtools view -@ ${task.cpus} -bS mapped_1.sam > mapped_1.bam
    rm -rf mapped_1.sam
    samtools view -@ ${task.cpus} -bS mapped_2.sam > mapped_2.bam
    rm -rf mapped_2.sam
    samtools fastq -@ ${task.cpus} mapped_1.bam > mapped_1.fastq
    rm -rf mapped_1.bam
    samtools fastq -@ ${task.cpus} mapped_2.bam > mapped_2.fastq
    rm -rf mapped_2.bam
    """
}


process run_optitype {
  label 'r5_4xlarge_500g'
  container  "${params.container.optitype}"
  publishDir "${params.outdir_run}/hla_type/", 
              pattern: "optitype_results/${sample_id}",
              mode: 'copy', 
              overwrite: true
  cpus 16
  memory '120 GB'
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
      } else { 
        if (params.bam_type == 'bam') {   // url + bam
          download_files_url_bam(get_sample_id.out.res_ch)
          fastq_out = bam_to_fastq(download_files_url_bam.out.res_ch)
        }
        else {  // url + cram
          download_files_url_cram(get_sample_id.out.res_ch)
          fastq_out = cram_to_fastq(download_files_url_cram.out.res_ch)
        }
      }
      reads_mapping(
        fastq_out.res_ch,
        params.hla_ref_prefix,
        Channel.fromPath(params.hla_ref).collect()
      )
      run_samtools(
        Channel.fromPath(params.hla_ref).collect(), 
        reads_mapping.out.res_ch
      )
     run_optitype(run_samtools.out.res_ch)

     emit:
       hla_typing_out = run_optitype.out.res_ch
}
