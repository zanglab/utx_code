#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WTK27M3.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab202012CHIP/20201209CHIP/"
outdir="./process_qc_out/WTK27M3"

bwa mem -t 8 ${genome_index} ${fqdir}/WTK27M3_S7_L001_R1_001.fastq.gz > ${outdir}/WTK27M3_WTK27M3_S7_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTK27M3_WTK27M3_S7_L001_R1_001.sam > ${outdir}/WTK27M3_WTK27M3_S7_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTK27M3_S7_L002_R1_001.fastq.gz > ${outdir}/WTK27M3_WTK27M3_S7_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTK27M3_WTK27M3_S7_L002_R1_001.sam > ${outdir}/WTK27M3_WTK27M3_S7_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTK27M3_S7_L003_R1_001.fastq.gz > ${outdir}/WTK27M3_WTK27M3_S7_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTK27M3_WTK27M3_S7_L003_R1_001.sam > ${outdir}/WTK27M3_WTK27M3_S7_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/WTK27M3_S7_L004_R1_001.fastq.gz > ${outdir}/WTK27M3_WTK27M3_S7_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WTK27M3_WTK27M3_S7_L004_R1_001.sam > ${outdir}/WTK27M3_WTK27M3_S7_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/WTK27M3_WTK27M3_S7_L001_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L002_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L003_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L004_R1_001.bam -n ${outdir}/WTK27M3

bdg2bw ${outdir}/WTK27M3_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WTK27M3_qc.log -s hg38 -b WTK27M3 -f .fastq.gz

samtools merge -@ 4 ${outdir}/WTK27M3_treat.bam \
${outdir}/WTK27M3_WTK27M3_S7_L001_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L002_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L003_R1_001.bam \
${outdir}/WTK27M3_WTK27M3_S7_L004_R1_001.bam

