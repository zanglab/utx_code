#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/MT2K27AC.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab-20201220seq"
outdir="./process_qc_out/MT2K27AC"

bwa mem -t 8 ${genome_index} ${fqdir}/MT2K27AC_S15_L001_R1_001.fastq.gz > ${outdir}/MT2K27AC_MT2K27AC_S15_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/MT2K27AC_MT2K27AC_S15_L001_R1_001.sam > ${outdir}/MT2K27AC_MT2K27AC_S15_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/MT2K27AC_S15_L002_R1_001.fastq.gz > ${outdir}/MT2K27AC_MT2K27AC_S15_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/MT2K27AC_MT2K27AC_S15_L002_R1_001.sam > ${outdir}/MT2K27AC_MT2K27AC_S15_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/MT2K27AC_S15_L003_R1_001.fastq.gz > ${outdir}/MT2K27AC_MT2K27AC_S15_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/MT2K27AC_MT2K27AC_S15_L003_R1_001.sam > ${outdir}/MT2K27AC_MT2K27AC_S15_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/MT2K27AC_S15_L004_R1_001.fastq.gz > ${outdir}/MT2K27AC_MT2K27AC_S15_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/MT2K27AC_MT2K27AC_S15_L004_R1_001.sam > ${outdir}/MT2K27AC_MT2K27AC_S15_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/MT2K27AC_MT2K27AC_S15_L001_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L002_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L003_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L004_R1_001.bam -n ${outdir}/MT2K27AC

bdg2bw ${outdir}/MT2K27AC_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/MT2K27AC_qc.log -s hg38 -b MT2K27AC -f .fastq.gz

samtools merge -@ 4 ${outdir}/MT2K27AC_treat.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L001_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L002_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L003_R1_001.bam \
${outdir}/MT2K27AC_MT2K27AC_S15_L004_R1_001.bam

