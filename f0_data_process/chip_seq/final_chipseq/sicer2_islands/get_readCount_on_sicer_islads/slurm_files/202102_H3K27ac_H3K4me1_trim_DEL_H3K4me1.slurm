#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_202102_H3K27ac_H3K4me1_trim_DEL_H3K4me1.out

python get_readCount_on_region_from_BamBed.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer2_islands/merged_islands/H3K4me1.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/DEL_H3K4me1/DEL_H3K4me1_treat.bam \
-s hg38 -f bam -o readCount_csv/202102_H3K27ac_H3K4me1_trim/DEL_H3K4me1_on_sicer_islands.csv

