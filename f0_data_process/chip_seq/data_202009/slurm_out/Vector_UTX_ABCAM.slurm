#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/Vector_UTX_ABCAM.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/BS202009_JiangLab/Fastq/"
outdir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202009/process_qc_out/Vector_UTX_ABCAM"

bwa mem -t 8 ${genome_index} ${fqdir}/BS9-PCDH-UTXABCAM_S5_L001_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L001_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS9-PCDH-UTXABCAM_S5_L002_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L002_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS9-PCDH-UTXABCAM_S5_L003_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L003_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS9-PCDH-UTXABCAM_S5_L004_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L004_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L004_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS17-PCDH-INPUT_S9_L001_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L001_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L001_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS17-PCDH-INPUT_S9_L002_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L002_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L002_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS17-PCDH-INPUT_S9_L003_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L003_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L003_R1_001.bam

bwa mem -t 8 ${genome_index} ${fqdir}/BS17-PCDH-INPUT_S9_L004_R1_001.fastq.gz > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L004_R1_001.sam > ${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L001_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L002_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L003_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L004_R1_001.bam -c \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L001_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L002_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L003_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L004_R1_001.bam -n ${outdir}/Vector_UTX_ABCAM

bdg2bw ${outdir}/Vector_UTX_ABCAM_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/Vector_UTX_ABCAM_qc.log -s hg38 -b Vector_UTX_ABCAM -f .fastq.gz

samtools merge -@ 4 ${outdir}/Vector_UTX_ABCAM_treat.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L001_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L002_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L003_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS9-PCDH-UTXABCAM_S5_L004_R1_001.bam

samtools merge -@ 4 ${outdir}/Vector_UTX_ABCAM_control.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L001_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L002_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L003_R1_001.bam \
${outdir}/Vector_UTX_ABCAM_BS17-PCDH-INPUT_S9_L004_R1_001.bam

