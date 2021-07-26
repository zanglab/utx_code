#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WTP300.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab-20201220seq"
outdir="./process_qc_out/WTP300"

bwa mem -t 8 ${genome_index} ${fqdir}/WTP300_S6_L001_R1_001.fastq.gz > ${outdir}/WTP300_WTP300_S6_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTP300_WTP300_S6_L001_R1_001.sam > ${outdir}/WTP300_WTP300_S6_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTP300_S6_L002_R1_001.fastq.gz > ${outdir}/WTP300_WTP300_S6_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTP300_WTP300_S6_L002_R1_001.sam > ${outdir}/WTP300_WTP300_S6_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTP300_S6_L003_R1_001.fastq.gz > ${outdir}/WTP300_WTP300_S6_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTP300_WTP300_S6_L003_R1_001.sam > ${outdir}/WTP300_WTP300_S6_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTP300_S6_L004_R1_001.fastq.gz > ${outdir}/WTP300_WTP300_S6_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTP300_WTP300_S6_L004_R1_001.sam > ${outdir}/WTP300_WTP300_S6_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/WTP300_WTP300_S6_L001_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L002_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L003_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L004_R1_001.bam -n ${outdir}/WTP300

bdg2bw ${outdir}/WTP300_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WTP300_qc.log -s hg38 -b WTP300 -f .fastq.gz

samtools merge -@ 4 ${outdir}/WTP300_treat.bam \
${outdir}/WTP300_WTP300_S6_L001_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L002_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L003_R1_001.bam \
${outdir}/WTP300_WTP300_S6_L004_R1_001.bam

