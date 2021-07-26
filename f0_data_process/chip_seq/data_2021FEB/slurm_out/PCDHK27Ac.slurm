#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/PCDHK27Ac.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHK27Ac_S10_L001_R1_001.fastq.gz > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L001_R1_001.sam > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHK27Ac_S10_L002_R1_001.fastq.gz > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L002_R1_001.sam > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHK27Ac_S10_L003_R1_001.fastq.gz > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L003_R1_001.sam > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021FEBCHIP//PCDHK27Ac_S10_L004_R1_001.fastq.gz > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L004_R1_001.sam > ./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L001_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L002_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L003_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L004_R1_001.bam -n ./process_qc_out/PCDHK27Ac/PCDHK27Ac

bdg2bw ./process_qc_out/PCDHK27Ac/PCDHK27Ac_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/PCDHK27Ac/PCDHK27Ac_treat.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L001_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L002_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L003_R1_001.bam \
./process_qc_out/PCDHK27Ac/PCDHK27Ac_PCDHK27Ac_S10_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/PCDHK27Ac_qc.log -s hg38 -b PCDHK27Ac -f .fastq.gz

