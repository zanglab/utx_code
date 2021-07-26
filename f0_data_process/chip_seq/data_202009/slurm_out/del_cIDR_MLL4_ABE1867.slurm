#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/del_cIDR_MLL4_ABE1867.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/BS202009_JiangLab/Fastq/"
outdir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202009/process_qc_out/del_cIDR_MLL4_ABE1867"

bwa mem -t 8 ${genome_index} ${fqdir}/BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS19-DEL3-INPUT_S11_L001_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L001_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS19-DEL3-INPUT_S11_L002_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L002_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS19-DEL3-INPUT_S11_L003_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L003_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS19-DEL3-INPUT_S11_L004_R1_001.fastq.gz > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L004_R1_001.sam > ${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.bam -c \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L001_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L002_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L003_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L004_R1_001.bam -n ${outdir}/del_cIDR_MLL4_ABE1867

bdg2bw ${outdir}/del_cIDR_MLL4_ABE1867_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/del_cIDR_MLL4_ABE1867_qc.log -s hg38 -b del_cIDR_MLL4_ABE1867 -f .fastq.gz

samtools merge -@ 4 ${outdir}/del_cIDR_MLL4_ABE1867_treat.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L001_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L002_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L003_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS23-DEL3-MLL4ABCAM_S15_L004_R1_001.bam

samtools merge -@ 4 ${outdir}/del_cIDR_MLL4_ABE1867_control.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L001_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L002_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L003_R1_001.bam \
${outdir}/del_cIDR_MLL4_ABE1867_BS19-DEL3-INPUT_S11_L004_R1_001.bam

