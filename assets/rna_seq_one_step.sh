#!/bin/sh

#SBATCH --job-name=rna_seq.sh
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -e -x

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
aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/reference/ucsc_refseq_hg38_20180629_genome/hg38 /tmp/reference/hg38 --recursive
aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/reference/ucsc_refseq_hg38_20180629_genome/genomeDir "/tmp/${sample_ID}/genomeDir" --recursive

mv "/tmp/${sample_ID}/${rna_r1_ID}/${rna_r1_name}" "/tmp/${sample_ID}/${sample_ID}_rna_R1.fastq.gz"
mv "/tmp/${sample_ID}/${rna_r2_ID}/${rna_r2_name}" "/tmp/${sample_ID}/${sample_ID}_rna_R2.fastq.gz"

rm -r "/tmp/${sample_ID}/${rna_r1_ID}/"
rm -r "/tmp/${sample_ID}/${rna_r2_ID}/"

docker run -u 510 -v "/tmp/${sample_ID}":/data/ zhanglab18/star:2.6.0a STAR --genomeDir /data/genomeDir \
                --runThreadN 8 \
                --outSAMtype BAM Unsorted \
                --readFilesCommand zcat \
                --readFilesIn "/data/${sample_ID}_rna_R1.fastq.gz" "/data/${sample_ID}_rna_R2.fastq.gz" \
                --outFileNamePrefix "/data/${sample_ID}."

docker run -u 510 -v /tmp/:/data/  zhanglab18/star:2.6.0a STAR --runMode genomeGenerate \
                --genomeDir "/data/${sample_ID}/genomeDir" \
                --genomeFastaFiles  /data/reference/hg38/ucsc_hg38.fa \
                --sjdbFileChrStartEnd "/data/${sample_ID}/${sample_ID}.SJ.out.tab" \
                --sjdbOverhang 75 \
                --runThreadN 8

docker run -u 510 -v /tmp/${sample_ID}/:/data/  zhanglab18/star:2.6.0a STAR --genomeDir /data/genomeDir \
                --runThreadN 8 \
                --outSAMtype BAM Unsorted \
                --readFilesCommand zcat \
                --readFilesIn "/data/${sample_ID}_rna_R1.fastq.gz" "/data/${sample_ID}_rna_R2.fastq.gz" \
                --outFileNamePrefix "/data/${sample_ID}.rna."

rm "/tmp/${sample_ID}/${sample_ID}.Aligned.out.bam" "/tmp/${sample_ID}/${sample_ID}.SJ.out.tab"
rm -r "/tmp/${sample_ID}/genomeDir"
rm "/tmp/${sample_ID}/${sample_ID}_rna_R1.fastq.gz" "/tmp/${sample_ID}/${sample_ID}_rna_R2.fastq.gz"

#aws s3 cp "/tmp/${sample_ID}" "s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/STAR/${sample_ID}" --recursive

#aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/STAR/${sample_ID} /tmp/${sample_ID} --recursive

docker run -u 510 -v "/tmp/${sample_ID}/":/data/ biocontainers/picard:2.3.0 java -Xmx16g -jar /opt/conda/share/picard-2.3.0-0/picard.jar AddOrReplaceReadGroups \
                I="/data/${sample_ID}.rna.Aligned.out.bam" \
                O="${sample_ID}.rna.Aligned.grouped.sorted.bam" \
                SO=coordinate \
                RGLB="${sample_ID}.rna-seq" \
                RGPL=illumina \
                RGPU="${sample_ID}.rna-seq" \
                RGSM="${sample_ID}.rna-seq"

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.out.bam"

#aws s3 cp /tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.bam s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/
#aws s3 cp /tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.output.metrics s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/

#aws s3 cp s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.bam /tmp/${sample_ID}/
docker run -u 510 -v "/tmp/${sample_ID}/":/data/ biocontainers/picard:2.3.0 java -Xmx16g -jar /opt/conda/share/picard-2.3.0-0/picard.jar MarkDuplicates \
                I="/data/${sample_ID}.rna.Aligned.grouped.sorted.bam" \
                O="/data/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.bam" \
                METRICS_FILE="/data/${sample_ID}.rna.Aligned.grouped.sorted.output.metrics" \
                VALIDATION_STRINGENCY=SILENT \
                CREATE_INDEX=true

docker run -u 510 -v "/tmp/${sample_ID}/":/data/ zhanglab18/samtools:1.9 samtools flagstat "/data/${sample_ID}.rna.Aligned.grouped.sorted.bam" > "${sample_ID}.rna.Aligned.grouped.sorted.mapping.statistic.txt"

mv "${sample_ID}.rna.Aligned.grouped.sorted.mapping.statistic.txt" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.mapping.statistic.txt"

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.bam"

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T SplitNCigarReads \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.bam" \
                -o "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.bam" \
                -RMQF 255 \
                -RMQT 60 \
                -rf ReassignOneMappingQuality \
                -U ALLOW_N_CIGAR_READS

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.bam" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.bai"

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T RealignerTargetCreator \
                -nt 8 \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.bam" \
                -known /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -known /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -o "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.intervals.list"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T IndelRealigner \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.bam" \
                -known /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -known /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -targetIntervals "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.intervals.list" \
                -o "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.bam"

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.bam" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.bai" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.intervals.list"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T BaseRecalibrator \
                -nct 8 \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.bam" \
                --knownSites /data/reference/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
                --knownSites /data/reference/hg38/dbsnp_138.hg38.vcf.gz \
                --knownSites /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                --knownSites /data/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -o "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.data.table"

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T PrintReads \
                -nct 8 \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.bam" \
                -BQSR "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.data.table" \
                -o "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.processed.bam"

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.bam" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.bai" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.data.table"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T HaplotypeCaller \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -I "/data/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.processed.bam" \
                --dbsnp /data/reference/hg38/dbsnp_138.hg38.vcf.gz \
                -dontUseSoftClippedBases \
                -stand_call_conf 20.0 \
                -o "/data/${sample_ID}/${sample_ID}.haplotypecaller.vcf"

rm "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.processed.bam" "/tmp/${sample_ID}/${sample_ID}.rna.Aligned.grouped.sorted.deduplicated.splitted.realigned.processed.bai"

docker run -u 510 -v /tmp/:/data/ biocontainers/picard:2.3.0 java -Xmx16g -jar /opt/conda/share/picard-2.3.0-0/picard.jar SortVcf \
                I= "/data/${sample_ID}/${sample_ID}.haplotypecaller.vcf" \
                O= "/data/${sample_ID}/${sample_ID}.haplotypecaller.sorted.vcf" \
                SEQUENCE_DICTIONARY=/data/reference/hg38/ucsc_hg38.dict

rm "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.vcf.idx" "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.vcf"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v /tmp/:/data/ broadinstitute/gatk3:3.8-0 java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R /data/reference/hg38/ucsc_hg38.fa \
                -V "/data/${sample_ID}/${sample_ID}.haplotypecaller.sorted.vcf" \
                -window 35 \
		-cluster 3 \
                -filterName FS \
                -filter "FS > 30.0" \
                -filterName QD \
                -filter "QD < 2.0" \
                -o "/data/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.vcf"

rm "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.vcf"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v "/tmp/${sample_ID}/":/data/ zhanglab18/vcftools:0.1.16 vcftools \
                --vcf "/data/${sample_ID}.haplotypecaller.sorted.filtered.vcf" \
                --remove-filtered-all \
                --recode \
                --out "/data/${sample_ID}.haplotypecaller.sorted.filtered"

rm "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.vcf"

perl /home/kail/software/annovar/table_annovar.pl "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.recode.vcf" /tmp/reference/hg38 \
                -buildver hg38 \
                -out "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.recode" \
                -protocol refGene \
                -operation g \
                -nastring . \
                --thread 8 \
                --maxgenethread 8 \
                -polish \
                -vcfinput

rm "/tmp/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.recode.vcf"

#aws s3 cp /tmp/${sample_ID}/ s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/${sample_ID}/ --recursive

docker run -u 510 -v /tmp/:/data/ zhanglab18/customprodbj:1.1.0 java -Xmx16g -jar /opt/customprodbj.jar \
                -i "/data/${sample_ID}/${sample_ID}.haplotypecaller.sorted.filtered.recode.hg38_multianno.txt" \
                -d /data/reference/hg38/hg38_refGeneMrna.fa \
                -r /data/reference/hg38/hg38_refGene.txt \
                -o "/data/${sample_ID}/" \
                -t

aws s3 cp "/tmp/${sample_ID}" "s3://zhanglab-kail/projects/2019_10_hnscc/data/results/rna_seq/final_results/${sample_ID}/" --recursive

rm -r "/tmp/${sample_ID}"
rm -r /tmp/reference

