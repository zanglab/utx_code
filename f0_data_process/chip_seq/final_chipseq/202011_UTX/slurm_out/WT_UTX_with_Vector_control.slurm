#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WT_UTX_with_Vector_control.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_L001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_L002.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_L003.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_L004.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTXNOV_S2_L001_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTXNOV_S2_L002_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTXNOV_S2_L003_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTXNOV_S2_L004_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_S9_L001_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_S9_L002_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_S9_L003_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/WTUTX_S9_L004_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTX_L001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTX_L002.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTX_L003.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTX_L004.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTXNOV_S1_L001_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTXNOV_S1_L002_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTXNOV_S1_L003_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX/PCDHUTXNOV_S1_L004_R1_001.fastq.gz > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.sam > ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.bam -c \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.bam -n ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.bam -c \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.bam -n ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_q01

bdg2bw ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_treat.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTXNOV_S2_L004_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_WTUTX_S9_L004_R1_001.bam

samtools merge -@ 4 ./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_control.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L002.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L003.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTX_L004.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L001_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L002_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L003_R1_001.bam \
./process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_PCDHUTXNOV_S1_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WT_UTX_with_Vector_control_qc.log -s hg38 -b WT_UTX_with_Vector_control -f .fastq.gz

