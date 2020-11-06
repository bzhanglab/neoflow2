
process pre_processing {
	container "${params.container.r_tidyverse}"

	input:
    path('manifest.txt')
    val start
    val end

	output:
    path "*-mapping_file.tsv", emit: res_ch

	script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    manifest <- read_tsv("manifest.txt")
    case_start <- ifelse($start == -1, 1, $start)
    case_end <- ifelse($end == -1, nrow(manifest), $end)
    stopifnot(case_start >= 1, case_start <= nrow(manifest))
    stopifnot(case_end >= 1, case_end <= nrow(manifest))
    stopifnot(case_start <= case_end)
    a <- manifest[case_start:case_end, ]
    experiment_names <- unique(a[,"experiment"] %>% pull(1))

    for(f in experiment_names){
      dat <- a %>% filter(experiment == f) 
      # for now the file type is "somatic"
      dat <- dat %>%  
             select(experiment, sample, vcf_file_path, vcf_file_name) %>%
             rename(file=vcf_file_name) %>% 
             add_column(file_type = "somatic")
      out_file <- paste(f,"-mapping_file.tsv",sep="")
      write_tsv(dat,out_file)
    }
    """
}

process get_exp_name {
  container "${params.container.ubuntu}"

  input:
     path(mapping_file)
  
  output:
     tuple env(exp_name), path(mapping_file), emit: res_ch
     
  """
  mf=${mapping_file}
  exp_name=\${mf%-mapping_file.tsv}
  """
}


process prepare_annovar {
   container "${params.container.ubuntu}"
   cpus 8

   input:
      path('annovar.tgz')
      path('annovar_anno.tgz')

    output:
      path 'annovar', emit: annovar_ch
      path 'annovar_anno', emit: annovar_anno_ch
    
    """
    tar xvfz annovar.tgz
    mkdir annovar_anno
    tar xvfz annovar_anno.tgz -C annovar_anno --strip-components 1
    """
}


process variant_annotation {
	container "${params.container.variant_annotation}"
  cpus 2
	publishDir "${params.outdir_run}/variant_annotation/",
              mode: 'copy',
              pattern: 'exp_*/*',
              overwrite: true

	input:
    tuple val(exp_name), path(mapping_file)
    path('vcf.tgz')
    path('annovar')
    path('annovar_anno')

	 output:
     tuple val(exp_name), path('*/*_anno.txt'), emit: anno_file
     path '*/*_multianno.txt', emit: v_multianno_file

	script:
    ofile = "${exp_name}_anno.txt"
    """
    mv ${mapping_file} mapping_file.tsv
    mkdir "exp_${exp_name}"
    tar xvfz vcf.tgz
    python /usr/local/bin/variant_annotation.py -i mapping_file.tsv \
      -d annovar_anno \
      -b ${params.annovar_buildver} \
      -c ${task.cpus} \
      -o ./ \
      -f ${ofile} \
      -a "annovar/table_annovar.pl"
    mv *_anno.txt "exp_${exp_name}"
    mv *_multianno.txt "exp_${exp_name}"
   """
}


process database_construct {
 	container "${params.container.neoflow}"
  cpus 4
  memory '12 GB'
	publishDir "${params.outdir_run}/customized_database/", 
             mode: 'copy', 
             pattern: 'exp_*/*',
             overwrite: true 

	input:
    tuple val(exp_name), path(anno_txt)
    path('annovar_anno')
    path('*')

	output:
    tuple val(exp_name), path('exp_*/*')
    tuple val(exp_name), path('*/*-var.fasta'), emit: var_ch
    tuple val(exp_name), path('*/*_anno-var.fasta'), emit: target_customized_db_fa 
    // this may not be necessary since ref.fasta is the same for all experiment
    tuple val(exp_name), path('ref.fasta'), emit: ref_ch
    tuple val(exp_name), path('*/experiment_varinfo/*.txt'), emit: exp_varinfo_ch 
    tuple val(exp_name), path('*/sample_varinfo/*.txt'), emit: sample_varinfo_ch 

	script:
    mrna_fa   = "annovar_anno/${params.annovar_buildver}_${params.annovar_protocol}Mrna.fa"
    gene_anno = "annovar_anno/${params.annovar_buildver}_${params.annovar_protocol}.txt"
    output_dir = "./exp_${exp_name}"
    """
    java -jar /opt/customprodbj.jar \
      -f ${anno_txt} \
      -d ${mrna_fa} \
       -r ${gene_anno} \
       -t \
       -o ${output_dir} \
       -p2 ${exp_name} \
       -ref ref.fasta
    cd ${output_dir}
    mkdir sample_varinfo
    mkdir experiment_varinfo
    mv ${exp_name}_*-varInfo.txt experiment_varinfo
    mv *-varInfo.txt sample_varinfo
    """
}

process split_exp_varinfo {
  container "${params.container.r_tidyverse}"
  cpus 1

  publishDir "${params.outdir_run}/customized_database/", 
             mode: 'copy', 
             pattern: 'exp_*/sample_varinfo/*-new-varInfo.txt',
             overwrite: true 
  input:
    path 'manifest.txt'
    tuple val(exp_name), path('*') 

  output:
    tuple val(exp_name),
          path('exp_*/sample_varinfo/*-new-varInfo.txt'),
          emit: new_sample_varinfo_ch 

  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  library(data.table)
  library(stringr)

  mani <- read_tsv("manifest.txt", col_types = cols(experiment = col_character()))
  samples <- mani %>% filter(experiment==${exp_name}) %>% pull(sample)

  merg_file <- paste0(${exp_name}, "_anno-varInfo.txt")
  merg_infor <- fread(merg_file)
  dir.create(paste0('exp_', ${exp_name}, "/sample_varinfo"), recursive=TRUE)
  for(sample in samples) {
    column_name <- as.name(paste0("sample:", sample))
    filter_infor <- merg_infor %>%
      filter(UQ(column_name) == 1)
    write.table(filter_infor, file = paste0("exp_", ${exp_name}, "/sample_varinfo/", sample, "-new-varInfo.txt"), 
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  """           
}

process format_db {
  container "${params.container.python}"
	publishDir "${params.outdir_run}/customized_database/",
             mode: 'copy',
             pattern: 'exp_*/*',
             overwrite: true
  cpus 8

	input:
    tuple val(exp_name),
          path('target_customized_db.fasta')

	output:
    tuple val(exp_name),
          path('*/*_format.fasta'),
          emit: out_ch

 	"""
  #!/usr/bin/env python
## https://github.com/bzhanglab/neoflow/blob/master/bin/format_db.py
import re
import sys
import os

in_db = 'target_customized_db.fasta'
out_db = re.sub(".fasta\$","_format.fasta",in_db)
print(str(out_db)+"\\n")
of = open(out_db,"w")
f = open(in_db,"r")
for line in f:
    if line.startswith(">"):
        of.write(line.split(" ")[0]+"\\n")
    else:
        of.write(line)
f.close()
of.close() 

exp_dir = 'exp_${exp_name}'
os.mkdir(exp_dir)
os.rename(out_db, os.path.join(exp_dir, out_db))
	"""
}


process generate_decoy_db{
  container "${params.container.pga}"
  cpus 1
	publishDir "${params.outdir_run}/customized_database/",
             mode: 'copy',
             pattern: 'exp_*/*',
             overwrite: true

	input:
    path('contaminants.fasta')
    tuple val(exp_name),
          path('format_db.fasta')

	output:
    tuple val(exp_name),
          path('*/*_target_decoy.fasta'), 
          emit: search_db_ch 

	script:
    """
    #!/usr/bin/env /usr/local/bin/Rscript
    library(PGA)
    output_dir <- file.path(paste0("exp_", ${exp_name}))
    dir.create(output_dir)
    out_target_decoy_db <- file.path(paste0("exp_", ${exp_name}),
         paste0(${exp_name}, "_target_decoy.fasta"))
    buildTargetDecoyDB(db="format_db.fasta",
      decoyPrefix="XXX_",
      cont_file="contaminants.fasta",
      output=out_target_decoy_db)

    """
}


workflow database_construction {
  main:
    pre_processing(params.manifest, params.start, params.end)
    prepare_annovar(
      params.annovar_file,
      params.annovar_anno_file
    )
    get_exp_name(pre_processing.out.res_ch.flatten())
    variant_annotation(
      get_exp_name.out.res_ch,
      params.vcf_file,
      prepare_annovar.out.annovar_ch,
      prepare_annovar.out.annovar_anno_ch
    )
    database_construct(
      variant_annotation.out.anno_file,
      prepare_annovar.out.annovar_anno_ch,
      variant_annotation.out.v_multianno_file
    )
    split_exp_varinfo(params.manifest, 
      database_construct.out.exp_varinfo_ch
    )
    format_db(
      database_construct.out.target_customized_db_fa
    )
    generate_decoy_db(
      params.contaminants,
      format_db.out.out_ch
    )
  
    emit:
      search_db_ch = generate_decoy_db.out.search_db_ch
      ref_ch = database_construct.out.ref_ch
      var_db_ch = database_construct.out.target_customized_db_fa
      sample_varinfo_ch = split_exp_varinfo.out.new_sample_varinfo_ch
      var_pep_info = database_construct.out.exp_varinfo_ch
}


