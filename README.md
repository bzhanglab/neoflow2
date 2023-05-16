# NeoFlow2: a proteogenomics pipeline for neoantigen discovery

NeoFlow2 is a streamlined computational workflow that integrates WES and MS/MS proteomics data for 
neoantigen prioritization to facilitate cancer immunotherapy. It includes four modules: 

  * Variant annotation and customized database construction
  * Variant peptide identification including MS/MS searching, FDR estimation, PepQuery validation
  * Human leukocyte antigen (HLA) typing
  * MHC-binding prediction and neoantigen prioritization. 

  The four modules are streamlined using [Nextflow](http://nextflow.io) and Docker. NeoFlow2 supports both label free and iTRAQ/TMT data. 
NeoFlow2 could be run in a single Linux computer or Could environment like AWS.

## How to run 
* run with docker locally

```console
nextflow run bzhanglab/neoflow2 -r main -profile docker \
   -params-file /path/to/my/parameter_file.yml
```

* run with aws batch on AWS cloud. For more information, please refer to
  [Running nextflow on AWS cloud](https://www.nextflow.io/docs/latest/awscloud.html).

```console
nextflow run bzhanglab/neoflow2 -r main -profile awsbatch \
   -bucket-dir s3://mybucket/workdir/2022-05-31 \
   -params-file /path/to/my/parameter_file.yml
```

To run with awsbatch, you must specify an s3 path as `outdir`, e.g.
`--outdir s3://mybucket/myfolder`.  In this case, results will be 
stored a directory under `s3://mybucket/myfolder`. The directory 
name is defined by the parameter `run_version`.


The parameter file is a yaml file that contains all the parameters. Below is a sample parameter file:

```
---
run_version: "2022_05_31_run"
outdir: "s3://mybucket/myfoler"
manifest: "/path/to/manifest/file.tsv"
maf: "/path/to/maf/file"
fusion_file: "/path/to/tared/gzipped/fusion/file"
bam_source: "uuid"
bam_type: "bam"
```

The following table lists the parameters and their default values if available.

| parameter name           | note                                                         |
| ------------------------ | ------------------------------------------------------------ |
| run_version              | name of the current run                                      |
| outdir                   | path to a local directory (/data/neoflow2_output) or s3 path (e.g. s3://mybucket/neoflow2_outpt) |
| manifest                 | path to input manifest file (see below for more detail)      |
| maf                      | path to maf file which contains somatic mutations for all the samples included in the manifest file |
| fusion_file              | path to tared and gzipped fusion file                        |
| bam_source               | `"uuid"` (gdc) or `"url"` (http) or `"path"` (s3, or local file path), Default is `"uuid"`|
| bam_type                 | `"bam"` or `"cram"`, Default is `"bam"`)                              |
| database                 | path to database output of previously run for global proteomics (this is used for running phosphoproteomics only) |
| hlatyping                | path to hla typing output of previously run for global proteomics (this is used for running phosphoproteomics only) |
| annovar_protocol         | The parameter of "protocol" for ANNOVAR, default is "refGene". Find more about the setting of this parameter at https://annovar.openbioinformatics.org/en/latest/user-guide/startup/ |
| annovar_anno_file        | ANNOVAR annotation datafile. All required annotation files are needed to be in a single TGZ format file |
| annovar_file             | ANNOVAR package file path. This is a TGZ format file which contains the ANNOVAR package. |
| annovar_buildver         | The genome build version. For example, hg19 or hg38. This is used for variant annotation using ANNOVAR for parameter "buildver". Find more about the setting of this parameter at https://annovar.openbioinformatics.org/en/latest/user-guide/startup/. Default is "GRCh38.p13.gencode.v34.basic". |
| hla_ref_prefix           | HLA DNA reference file in FASTA format. Default is "hla_reference_dna.fasta" |
| hla_ref                  | HLA DNA reference file path, e.g.: "/path/to/hla_reference/hla_reference_dna.fasta*". |
| seqtype                  | Reads type, "dna" or "rna". Default is "dna".                |
| search_engine            | The search engine used for MS/MS searching, comet=Comet, msgf=MS-GF+ or xtandem=X!Tandem. Default is "msgf". |
| search_para_file         | Parameter file for the search engine used in MS/MS searching. For MS-GF+, how to set the parameter file could be find at https://msgfplus.github.io/msgfplus/MSGFPlus.html and this parameter file is used by the parameter "-conf" in MS-GF+. For Comet, how to set the parameter file could be find at http://comet-ms.sourceforge.net/parameters/ and this parameter file is used by the parameter "-P" in Comet. For X!Tandem, the parameter file setting could be found at https://www.thegpm.org/TANDEM/api/index.html |
| search_engine_mem        | The memory limitation for MS/MS searching. Default is 60. The unit is GB. |
| contaminants             | Contaminant protein sequence file in FASTA format. This file could be downloaded from MaxQuant website |
| pv_enzyme                | Enzyme used for protein digestion. 0:Non enzyme, 1:Trypsin (default), 2:Trypsin (no P rule), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C.This is used by PepQuery. |
| pv_c                     | The max missed cleavages, default is 2. This is used by PepQuery. |
| pv_tol                   | Precursor ion m/z tolerance, default is 10. This is used by PepQuery. |
| pv_tolu                  | The unit of --tol, ppm or Da. Default is ppm. This is used by PepQuery. |
| pv_itol                  | The error window for fragment ion, default is 0.05. This is used by PepQuery. |
| pv_fixmod                | Fixed modification. A list of numbers separated by commas. E.g."1,2,3". Different modification is represented by different number. This is used by PepQuery. Default: "108,89,6" |
| pv_varmod                | Variable modification. The format is the same as pv_fixmod. This is used by PepQuery. Default: "117" |
| netmhc_file              | NetMHCpan 4.0 package in TGZ format.                         |
| neoantigen_output_prefix | The output folder name. Default is "neoflow".                |
| pga_prefix               | The prefix of PGA output file. Default is "pga".             |

## Input manifest file format (tsv)

The input manifest file should contain the following columns:
* sample
* experiment
* wxs_file_name   
* wxs_file_location
* mzml_files
  * a list of mzml files for each experiment, separated by comma
  * samples in the same experiment has the same `mzml_files` value
* mzml_path
  * path to mzml file for each experiment, provided as a tar file
  * samples in the same experiment has the same `mzml_path` value
* fusion (if fusion information is not available, 'NA' should be used)
  * path to the fusion tsv file for each sample
  * the path should match the path in the gzipped tar file (see next item)
  * all the fusion files should be tared and gzipped into a single file and provided as 
    a command line parameter (`--fusion_file`)

Here is a [sample manifest file](https://github.com/bzhanglab/neoflow2/blob/main/sample_manifest.txt).


## Output

The output of NeoFlow2 is organized into folders shown below:

<img src="https://github.com/bzhanglab/neoflow2/blob/main/neoflow2_output.png" alt="neoflow2 output"/>

Below is an overview of the data in each folder:
* variant_annotation: Variant annotation result for each sample
* customized_database: Customized protein database for each sample (Label free data) or each TMT/iTRAQ experiment.
* msms_searching: MS/MS searching result. This folder contains the peptide identification files directly generated by a search engine.
* fdr_estimation: FDR estimation for peptide identification from a search engine
* pepquery: Novel peptide validation result using PepQuery for each sample.
* novel_peptide_identification: Identified novel peptides passed PepQuery validation.
* hla_type: HLA typing result for each sample.
* binding_prediction: MHC-peptide binding prediction result for each sample.
* neoantigen_prediction: This folder contains the final output files of neoantigen prediction.

For each sample, the final output file is a text file in tab delimited format. Below is the description of the columns in the file:

| Columns                       | Description                                  |
|-------------------------------|----------------------------------------------|
| Variant_ID                    | variant ID defined by neoflow                |
| Chr                           | variant chromosome                           |
| Start                         | start position on genome                     |
| End                           | end position on genome                       |
| Ref                           | reference base                               |
| Alt                           | alterative base                              |
| Variant_Type                  | variant type annotated by ANNOVAR            |
| Variant_Function              | variant function annotated by ANNOVAR        |
| Gene                          | gene ID                                      |
| mRNA                          | mRNA ID                                      |
| Neoepitope                    | neoepitope peptide                           |
| Variant_Start                 | variant start position on neoepitope peptide |
| Variant_End                   | variant end position on neoepitope peptide   |
| AA_before                     | reference amino acid                         |
| AA_after                      | alterative amino acid                        |
| HLA_type                      | HLA type                                     |
| netMHCpan_binding_affinity_nM | MHC-peptide binding affinity                 |
| netMHCpan_precentail_rank     | MHC-peptide binding affinity rank            |
| protein_var_evidence_pep      | variant peptide                              |