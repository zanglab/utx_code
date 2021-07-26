#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/CRPCDHUTX.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab202012CHIP/202012CUTRUN/"
outdir="./process_qc_out/CRPCDHUTX"

bwa mem -t 8 ${genome_index} ${fqdir}/CRPCDHUTX_S1_L001_R1_001.fastq.gz > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L001_R1_001.sam > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRPCDHUTX_S1_L002_R1_001.fastq.gz > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L002_R1_001.sam > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRPCDHUTX_S1_L003_R1_001.fastq.gz > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L003_R1_001.sam > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRPCDHUTX_S1_L004_R1_001.fastq.gz > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L004_R1_001.sam > ${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L001_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L002_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L003_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L004_R1_001.bam -n ${outdir}/CRPCDHUTX

bdg2bw ${outdir}/CRPCDHUTX_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/CRPCDHUTX_qc.log -s hg38 -b CRPCDHUTX -f .fastq.gz

samtools merge -@ 4 ${outdir}/CRPCDHUTX_treat.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L001_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L002_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L003_R1_001.bam \
${outdir}/CRPCDHUTX_CRPCDHUTX_S1_L004_R1_001.bam

