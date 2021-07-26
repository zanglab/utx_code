#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/UTX_eIFIDR_MLL4_ABE1867.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/BS202009_JiangLab/Fastq/"
outdir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202009/process_qc_out/UTX_eIFIDR_MLL4_ABE1867"

bwa mem -t 8 ${genome_index} ${fqdir}/BS24-EIF-MLL4ABCAM_S16_L001_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L001_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS24-EIF-MLL4ABCAM_S16_L002_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L002_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS24-EIF-MLL4ABCAM_S16_L003_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L003_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS24-EIF-MLL4ABCAM_S16_L004_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L004_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L004_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS20-EIF-INPUT_S12_L001_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L001_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS20-EIF-INPUT_S12_L002_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L002_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS20-EIF-INPUT_S12_L003_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L003_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS20-EIF-INPUT_S12_L004_R1_001.fastq.gz > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L004_R1_001.sam > ${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L001_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L002_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L003_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L004_R1_001.bam -c \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L001_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L002_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L003_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L004_R1_001.bam -n ${outdir}/UTX_eIFIDR_MLL4_ABE1867

bdg2bw ${outdir}/UTX_eIFIDR_MLL4_ABE1867_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/UTX_eIFIDR_MLL4_ABE1867_qc.log -s hg38 -b UTX_eIFIDR_MLL4_ABE1867 -f .fastq.gz

samtools merge -@ 4 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_treat.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L001_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L002_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L003_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS24-EIF-MLL4ABCAM_S16_L004_R1_001.bam

samtools merge -@ 4 ${outdir}/UTX_eIFIDR_MLL4_ABE1867_control.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L001_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L002_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L003_R1_001.bam \
${outdir}/UTX_eIFIDR_MLL4_ABE1867_BS20-EIF-INPUT_S12_L004_R1_001.bam

