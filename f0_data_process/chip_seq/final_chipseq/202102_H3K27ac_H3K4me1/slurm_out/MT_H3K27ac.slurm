#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/MT_H3K27ac.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27Ac_S14_L001_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L001_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27Ac_S14_L002_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L002_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27Ac_S14_L003_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L003_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27Ac_S14_L004_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L004_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L004_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27AC_S18_L001_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L001_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27AC_S18_L002_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L002_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27AC_S18_L003_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L003_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1/MT2K27AC_S18_L004_R1_001.fastq.gz > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L004_R1_001.sam > ./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L001_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L002_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L003_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L004_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L001_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L002_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L003_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L004_R1_001.bam -n ./process_qc_out/MT_H3K27ac/MT_H3K27ac

bdg2bw ./process_qc_out/MT_H3K27ac/MT_H3K27ac_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/MT_H3K27ac/MT_H3K27ac_treat.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L001_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L002_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L003_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27Ac_S14_L004_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L001_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L002_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L003_R1_001.bam \
./process_qc_out/MT_H3K27ac/MT_H3K27ac_MT2K27AC_S18_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/MT_H3K27ac_qc.log -s hg38 -b MT_H3K27ac -f .fastq.gz

