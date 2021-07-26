#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WT_UTX.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_L001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_L002_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_L003_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_L004_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTXNOV_S2_L001_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTXNOV_S2_L002_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTXNOV_S2_L003_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTXNOV_S2_L004_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_S9_L001_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_S9_L002_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_S9_L003_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim/202011_UTX_trim/trim_out/WTUTX_S9_L004_R1_001_trimmed.fq.gz > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.sam > ./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.bam -n ./process_qc_out/WT_UTX/WT_UTX

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.bam -n ./process_qc_out/WT_UTX/WT_UTX_q01

bdg2bw ./process_qc_out/WT_UTX/WT_UTX_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/WT_UTX/WT_UTX_treat.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L002_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L003_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_L004_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTXNOV_S2_L004_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L001_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L002_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L003_R1_001_trimmed.bam \
./process_qc_out/WT_UTX/WT_UTX_WTUTX_S9_L004_R1_001_trimmed.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/WT_UTX_qc.log -s hg38 -b WT_UTX -f .fq.gz

