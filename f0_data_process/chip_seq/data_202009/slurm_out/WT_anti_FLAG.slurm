#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WT_anti_FLAG.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/BS202009_JiangLab/Fastq/"
outdir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202009/process_qc_out/WT_anti_FLAG"

bwa mem -t 8 ${genome_index} ${fqdir}/BS2-FL-M2_S2_L001_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L001_R1_001.sam > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS2-FL-M2_S2_L002_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L002_R1_001.sam > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS2-FL-M2_S2_L003_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L003_R1_001.sam > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS2-FL-M2_S2_L004_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L004_R1_001.sam > ${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L004_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS18-FL-INPUT_S10_L001_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L001_R1_001.sam > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS18-FL-INPUT_S10_L002_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L002_R1_001.sam > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS18-FL-INPUT_S10_L003_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L003_R1_001.sam > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS18-FL-INPUT_S10_L004_R1_001.fastq.gz > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L004_R1_001.sam > ${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L001_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L002_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L003_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L004_R1_001.bam -c \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L001_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L002_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L003_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L004_R1_001.bam -n ${outdir}/WT_anti_FLAG

bdg2bw ${outdir}/WT_anti_FLAG_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WT_anti_FLAG_qc.log -s hg38 -b WT_anti_FLAG -f .fastq.gz

samtools merge -@ 4 ${outdir}/WT_anti_FLAG_treat.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L001_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L002_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L003_R1_001.bam \
${outdir}/WT_anti_FLAG_BS2-FL-M2_S2_L004_R1_001.bam

samtools merge -@ 4 ${outdir}/WT_anti_FLAG_control.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L001_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L002_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L003_R1_001.bam \
${outdir}/WT_anti_FLAG_BS18-FL-INPUT_S10_L004_R1_001.bam

