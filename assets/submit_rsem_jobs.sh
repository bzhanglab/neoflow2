#!/bin/sh

set -e -x

table=~/projects/2019_10_hnscc/conf/data_tables/hnscc_tumor_sample_table.txt

sed 1d $table | while read -r Proteomics_Participant_ID tumor_bam tumor_bam_id tumor_bam_md5 RNAseq_R1_filename RNAseq_R2_filename RNAseq_R1_UUID RNAseq_R2_UUID
do
	echo $Proteomics_Participant_ID
	sbatch rna_rsem.sh $Proteomics_Participant_ID $RNAseq_R1_UUID $RNAseq_R1_filename $RNAseq_R2_UUID $RNAseq_R2_filename
done
