#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_re_1st_submission_H3K4me3_MLL4SC_trim_EIF_MLL4.out

python get_readCount_on_region_from_BamBed.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer2_islands/merged_islands/MLL4.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out/EIF_MLL4/EIF_MLL4_treat.bam \
-s hg38 -f bam -o readCount_csv/re_1st_submission_H3K4me3_MLL4SC_trim/EIF_MLL4_on_sicer_islands.csv

