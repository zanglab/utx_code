#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RP_slurms/slurm_Vector_H3K4me2_on_Promoter.out

python get_RP_from_readcount.py \
-i /Volumes/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_202012_H3K4me2_trim/process_qc_out/Vector_H3K4me2/Vector_H3K4me2_treat.bam \
-s hg38 -f bam -o RP_csv/Vector_H3K4me2_on_Promoter.csv

