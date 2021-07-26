#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/EIFK4M1.out 

bwa aln  -q 5 -l 32 -k 2 -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L001_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.sai
bwa samse /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.sai /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L001_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.sam > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.bam

bwa aln  -q 5 -l 32 -k 2 -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L002_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.sai
bwa samse /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.sai /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L002_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.sam > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.bam

bwa aln  -q 5 -l 32 -k 2 -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L003_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.sai
bwa samse /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.sai /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L003_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.sam > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.bam

bwa aln  -q 5 -l 32 -k 2 -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L004_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.sai
bwa samse /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.sai /nv/vol190/zanglab/bs6yk/202101CHIPSEQ//EIFK4M1_S16_L004_R1_001.fastq.gz > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.sam > ./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.bam -n ./process_qc_out/EIFK4M1/EIFK4M1

bdg2bw ./process_qc_out/EIFK4M1/EIFK4M1_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/EIFK4M1/EIFK4M1_treat.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L001_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L002_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L003_R1_001.bam \
./process_qc_out/EIFK4M1/EIFK4M1_EIFK4M1_S16_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/EIFK4M1_qc.log -s hg38 -b EIFK4M1 -f .fastq.gz

