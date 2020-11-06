NeoFlow2: a proteogenomics pipeline for neoantigen discovery


* run with docker locally

```console
nextflow run bzhanglab/neoflow2 -r main -profile docker \
   --manifest /path/to/my_manifest.tsv \
   --run_version neoflow_2020_11_06
```

* run with aws batch on AWS cloud

```console
nextflow run bzhanglab/neoflow2 -r main -profile awsbatch \
   --manifest /path/to/my_manifest.tsv \
   -bucket-dir s3://mybucket/workdir/2020-11-06 \
   --outdir s3://mybucket/neoflow2-results \
   --run_version neoflow_2020_11_06
```

To run with awsbatch, you must specify an s3 path as `outdir`, e.g.
`--outdir s3://mybucket/myfolder`.  In this case, results will be 
stored under `s3://mybucket/myfolder/run_version`.


Input manifest file format (tsv)

```
sample	experiment	wxs_file_name	wxs_file_uuid	mgf_file_name	mgf_file_path	vcf_file_name	vcf_file_path
C3L-00997	1	c123a154-adec-48b6-923a-9bc5a29d1601_wxs_gdc_realn.bam	38dfabd0-45ce-4434-80af-0d625e9a0cba	01CPTAC_HNSCC_Proteome_JHU_20190507.mgf	s3://zhanglab-kail/projects/2019_10_hnscc/data/mgf_files/01CPTAC_HNSCC_Proteome_JHU_20190507.mgf	C3L-00997_T.avinput	s3://zhanglab-kail/projects/neoflow2_development/vcf_files/C3L-00997_T.avinput
C3N-03849	1	58d317f1-0935-498a-8d5f-b16499308ab4_wxs_gdc_realn.bam	130912fb-80ca-42cd-bd76-284e3c843ea1	01CPTAC_HNSCC_Proteome_JHU_20190507.mgf	s3://zhanglab-kail/projects/2019_10_hnscc/data/mgf_files/01CPTAC_HNSCC_Proteome_JHU_20190507.mgf	C3N-03849_T.avinput	s3://zhanglab-kail/projects/neoflow2_development/vcf_files/C3N-03849_T.avinput
C3N-03664	2	3908eb72-303d-4316-b8a6-41f76e1bbed5_wxs_gdc_realn.bam	ec399f46-c282-4923-8744-252b8bc15d3c	02CPTAC_HNSCC_Proteome_JHU_20190513.mgf	s3://zhanglab-kail/projects/2019_10_hnscc/data/mgf_files/02CPTAC_HNSCC_Proteome_JHU_20190513.mgf	C3N-03664_T.avinput	s3://zhanglab-kail/projects/neoflow2_development/vcf_files/C3N-03664_T.avinput
C3N-03781	2	4fcb737a-a3a1-47b5-b8c7-9363fd57fd4d_wxs_gdc_realn.bam	a71b51a8-06f9-45a8-a7cb-3cbedbf72b66	02CPTAC_HNSCC_Proteome_JHU_20190513.mgf	s3://zhanglab-kail/projects/2019_10_hnscc/data/mgf_files/02CPTAC_HNSCC_Proteome_JHU_20190513.mgf	C3N-03781_T.avinput	s3://zhanglab-kail/projects/neoflow2_development/vcf_files/C3N-03781_T.avinpu
```

Important parameters:

```
  hla_ref_prefix: default 'hla_reference_dna.fasta'
  hla_ref: path to hla reference file
  annovar_buildver: default 'hg38'
  annovar_protocol: default'refGene'
  annovar_anno_file: path to annovar annotation file (gzipped  tar)
  seqtype: default "dna"
  annovar_file: path to annovar software (gzipped tar)
  vcf_file: path to all vcf files (gzipped tar)
  contaminants: path to the contamiants file (fasta)
  pv_enzyme: default 1
  pv_c: default  2
  pv_tol: default 10
  pv_tolu: default 'ppm'
  pv_itol: default 0.5
  pv_fixmod: default 6
  pv_varmod: default 107
  pga_prefix: default 'pga'
  // chose one of the search engine and set the path to its parameter file
  search_engine: 'comet'
  search_para_file: path to comet parameter file
  // search_engine: 'msgf'
  // search_para_file: path to msgf parameter file
  // search_engine:  'xtandem'
  // search_para_file: path to xtandem parameter file
  netmhc_file: path to metmhcpan software (gzipped tar)
  neoantigen_output_prefix: default 'neoflow'
```


Most of them can be set as the default values provided. Others can be 
provided through the command line with the format:

`--PARAMETER_NAME PARAMETER_VALUE`

