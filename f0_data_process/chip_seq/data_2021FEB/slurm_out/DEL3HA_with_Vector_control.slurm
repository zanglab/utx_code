#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/DEL3HA_with_Vector_control.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//DEL3HA_S3_L001_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L001_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//DEL3HA_S3_L002_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L002_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//DEL3HA_S3_L003_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L003_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//DEL3HA_S3_L004_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L004_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L004_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHHA_S1_L001_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L001_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHHA_S1_L002_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L002_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHHA_S1_L003_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L003_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHHA_S1_L004_R1_001.fastq.gz > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L004_R1_001.sam > ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L001_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L002_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L003_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L004_R1_001.bam -c \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L001_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L002_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L003_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L004_R1_001.bam -n ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control

bdg2bw ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_treat.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L001_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L002_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L003_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_DEL3HA_S3_L004_R1_001.bam

samtools merge -@ 4 ./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_control.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L001_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L002_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L003_R1_001.bam \
./process_qc_out/DEL3HA_with_Vector_control/DEL3HA_with_Vector_control_PCDHHA_S1_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/DEL3HA_with_Vector_control_qc.log -s hg38 -b DEL3HA_with_Vector_control -f .fastq.gz

