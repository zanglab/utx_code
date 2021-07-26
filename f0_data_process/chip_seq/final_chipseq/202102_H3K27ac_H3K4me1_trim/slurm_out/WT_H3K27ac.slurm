#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WT_H3K27ac.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/FLK27Ac_S11_L001_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L001_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L001_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L001_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/FLK27Ac_S11_L002_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L002_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L002_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L002_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/FLK27Ac_S11_L003_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L003_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L003_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L003_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/FLK27Ac_S11_L004_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L004_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L004_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L004_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/WTK27AC_S15_L001_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L001_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L001_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L001_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/WTK27AC_S15_L002_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L002_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L002_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L002_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/WTK27AC_S15_L003_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L003_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L003_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L003_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202102_H3K27ac_H3K4me1_trim/trim_out/WTK27AC_S15_L004_R1_001_trimmed.fq.gz > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L004_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L004_R1_001_trimmed.sam > ./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L004_R1_001_trimmed.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L001_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L002_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L003_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L004_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L001_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L002_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L003_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L004_R1_001_trimmed.bam -n ./process_qc_out/WT_H3K27ac/WT_H3K27ac

bdg2bw ./process_qc_out/WT_H3K27ac/WT_H3K27ac_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/WT_H3K27ac/WT_H3K27ac_treat.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L001_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L002_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L003_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_FLK27Ac_S11_L004_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L001_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L002_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L003_R1_001_trimmed.bam \
./process_qc_out/WT_H3K27ac/WT_H3K27ac_WTK27AC_S15_L004_R1_001_trimmed.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WT_H3K27ac_qc.log -s hg38 -b WT_H3K27ac -f .fq.gz

