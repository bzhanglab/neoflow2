#!/bin/sh

#SBATCH --job-name=rna_rsem
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -e -x

date

alias gdcz="/home/kail/software/gdc/gdc-client download -t /home/kail/software/gdc/gdc-user-token.2019-10-01T16_22_05.675Z.txt -d "
sample_ID=${1}
rna_r1_ID=${2}
rna_r1_name=${3}
rna_r2_ID=${4}
rna_r2_name=${5}

if [ -d "/tmp/${sample_ID}" ];then
        rm -r "/tmp/${sample_ID}"
        mkdir "/tmp/${sample_ID}"
else
        mkdir "/tmp/${sample_ID}"
fi

if [ -d /tmp/reference ];then
        rm -r /tmp/reference
        mkdir /tmp/reference
else
        mkdir /tmp/reference
fi

gdcz "/tmp/${sample_ID}" ${rna_r1_ID}
gdcz "/tmp/${sample_ID}" ${rna_r2_ID}

mv "/tmp/${sample_ID}/${rna_r1_ID}/${rna_r1_name}" "/tmp/${sample_ID}/${sample_ID}_rna_R1.fastq.gz"
mv "/tmp/${sample_ID}/${rna_r2_ID}/${rna_r2_name}" "/tmp/${sample_ID}/${sample_ID}_rna_R2.fastq.gz"

rm -r "/tmp/${sample_ID}/${rna_r1_ID}/"
rm -r "/tmp/${sample_ID}/${rna_r2_ID}/"

mkdir /tmp/${sample_ID}/genomeDir

aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/reference/ucsc_refseq_hg38_20180629_genome/hg38/ /tmp/reference --recursive
aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/reference/rna_seq/genomeDir /tmp/${sample_ID}/genomeDir --recursive

#docker run -u 510 -v /tmp/:/data/ zhanglab18/rsem:1.2.26 /opt/RSEM-1.2.26/rsem-prepare-reference \
#	--gtf /data/reference/ucsc_refseq_hg38_20180629_with_gene_symbol_for_rsem.gtf \
#	--star \
#	--star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
#	-p 8 \
#	/data/reference/ucsc_hg38.fa \
#	/data/${sample_ID}/genomeDir/rsem

#aws s3 cp /tmp/${sample_ID}/genomeDir/ s3://zhanglab-kail/projects/2019_10_hnscc/reference/rna_seq/genomeDir_new/ --recursive

docker run -u 510 -v /tmp/:/data/ zhanglab18/rsem:1.2.26 /opt/RSEM-1.2.26/rsem-calculate-expression \
                        --paired-end \
                        --star \
                        --star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
                        -p 8 \
                        --gzipped-read-file \
                        "/data/${sample_ID}/${sample_ID}_rna_R1.fastq.gz" \
                        "/data/${sample_ID}/${sample_ID}_rna_R2.fastq.gz" \
                        "/data/${sample_ID}/genomeDir/rsem" \
                        "/data/${sample_ID}/${sample_ID}"

rm "/tmp/${sample_ID}/${sample_ID}_rna_R1.fastq.gz" "/tmp/${sample_ID}/${sample_ID}_rna_R2.fastq.gz"
rm "/tmp/${sample_ID}/${sample_ID}.transcript.sorted.bam" "/tmp/${sample_ID}/${sample_ID}.transcript.sorted.bam.bai"
rm "/tmp/${sample_ID}/${sample_ID}.transcript.bam"
rm -r "/tmp/${sample_ID}/genomeDir"

aws s3 cp "/tmp/${sample_ID}" s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/rsem_results/${sample_ID} --recursive

date

rm -r "/tmp/${sample_ID}"
rm -r /tmp/reference
