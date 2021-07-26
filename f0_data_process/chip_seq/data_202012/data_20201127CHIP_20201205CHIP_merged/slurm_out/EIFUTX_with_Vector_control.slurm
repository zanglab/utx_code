#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/EIFUTX_with_Vector_control.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/chipseq/chipseq202012/"
outdir="./process_qc_out/EIFUTX_with_Vector_control"

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_L001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L001.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_L002.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L002.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L002.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L002.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_L003.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L003.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L003.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L003.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_L004.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L004.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L004.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_L004.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_S11_L001_R1_001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L001_R1_001.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_S11_L002_R1_001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L002_R1_001.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_S11_L003_R1_001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L003_R1_001.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/EIFUTX_S11_L004_R1_001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L004_R1_001.sam > ${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L004_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L001.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L001.sam > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L002.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L002.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L002.sam > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L002.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L003.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L003.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L003.sam > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L003.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L004.fastq.gz > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L004.sam
samtools view -bS -q 1 -@ 8 ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L004.sam > ${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L004.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L002.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L003.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L004.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L001_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L002_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L003_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L004_R1_001.bam -c \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L001.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L002.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L003.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L004.bam -n ${outdir}/EIFUTX_with_Vector_control

bdg2bw ${outdir}/EIFUTX_with_Vector_control_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/EIFUTX_with_Vector_control_qc.log -s hg38 -b EIFUTX_with_Vector_control -f .fastq.gz

samtools merge -@ 4 ${outdir}/EIFUTX_with_Vector_control_treat.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L002.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L003.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_L004.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L001_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L002_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L003_R1_001.bam \
${outdir}/EIFUTX_with_Vector_control_EIFUTX_S11_L004_R1_001.bam

samtools merge -@ 4 ${outdir}/EIFUTX_with_Vector_control_control.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L001.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L002.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L003.bam \
${outdir}/EIFUTX_with_Vector_control_PCDHUTX_L004.bam

