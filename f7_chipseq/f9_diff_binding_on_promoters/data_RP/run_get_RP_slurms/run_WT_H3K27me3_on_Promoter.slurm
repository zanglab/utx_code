#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RP_slurms/slurm_WT_H3K27me3_on_Promoter.out

python get_RP_from_readcount.py \
-i /Volumes/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/WT_H3K27me3/WT_H3K27me3_treat.bam \
-s hg38 -f bam -o RP_csv/WT_H3K27me3_on_Promoter.csv

