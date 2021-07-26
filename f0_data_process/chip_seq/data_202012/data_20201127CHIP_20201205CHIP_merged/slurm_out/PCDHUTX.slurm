#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/PCDHUTX.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/chipseq/chipseq202012/"
outdir="./process_qc_out/PCDHUTX"

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L001.fastq.gz > ${outdir}/PCDHUTX_PCDHUTX_L001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/PCDHUTX_PCDHUTX_L001.sam > ${outdir}/PCDHUTX_PCDHUTX_L001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L002.fastq.gz > ${outdir}/PCDHUTX_PCDHUTX_L002.sam
samtools view -bS -q 1 -@ 8 ${outdir}/PCDHUTX_PCDHUTX_L002.sam > ${outdir}/PCDHUTX_PCDHUTX_L002.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L003.fastq.gz > ${outdir}/PCDHUTX_PCDHUTX_L003.sam
samtools view -bS -q 1 -@ 8 ${outdir}/PCDHUTX_PCDHUTX_L003.sam > ${outdir}/PCDHUTX_PCDHUTX_L003.bam

bwa mem -t 8 ${genome_index} ${fqdir}/PCDHUTX_L004.fastq.gz > ${outdir}/PCDHUTX_PCDHUTX_L004.sam
samtools view -bS -q 1 -@ 8 ${outdir}/PCDHUTX_PCDHUTX_L004.sam > ${outdir}/PCDHUTX_PCDHUTX_L004.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/PCDHUTX_PCDHUTX_L001.bam \
${outdir}/PCDHUTX_PCDHUTX_L002.bam \
${outdir}/PCDHUTX_PCDHUTX_L003.bam \
${outdir}/PCDHUTX_PCDHUTX_L004.bam -n ${outdir}/PCDHUTX

bdg2bw ${outdir}/PCDHUTX_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/PCDHUTX_qc.log -s hg38 -b PCDHUTX -f .fastq.gz

samtools merge -@ 4 ${outdir}/PCDHUTX_treat.bam \
${outdir}/PCDHUTX_PCDHUTX_L001.bam \
${outdir}/PCDHUTX_PCDHUTX_L002.bam \
${outdir}/PCDHUTX_PCDHUTX_L003.bam \
${outdir}/PCDHUTX_PCDHUTX_L004.bam

