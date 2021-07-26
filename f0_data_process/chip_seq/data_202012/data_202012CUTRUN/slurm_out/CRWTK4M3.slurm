#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/CRWTK4M3.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab202012CHIP/202012CUTRUN/"
outdir="./process_qc_out/CRWTK4M3"

bwa mem -t 8 ${genome_index} ${fqdir}/CRWTK4M3_S7_L001_R1_001.fastq.gz > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRWTK4M3_CRWTK4M3_S7_L001_R1_001.sam > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRWTK4M3_S7_L002_R1_001.fastq.gz > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRWTK4M3_CRWTK4M3_S7_L002_R1_001.sam > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRWTK4M3_S7_L003_R1_001.fastq.gz > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRWTK4M3_CRWTK4M3_S7_L003_R1_001.sam > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/CRWTK4M3_S7_L004_R1_001.fastq.gz > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/CRWTK4M3_CRWTK4M3_S7_L004_R1_001.sam > ${outdir}/CRWTK4M3_CRWTK4M3_S7_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L001_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L002_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L003_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L004_R1_001.bam -n ${outdir}/CRWTK4M3

bdg2bw ${outdir}/CRWTK4M3_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/CRWTK4M3_qc.log -s hg38 -b CRWTK4M3 -f .fastq.gz

samtools merge -@ 4 ${outdir}/CRWTK4M3_treat.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L001_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L002_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L003_R1_001.bam \
${outdir}/CRWTK4M3_CRWTK4M3_S7_L004_R1_001.bam

