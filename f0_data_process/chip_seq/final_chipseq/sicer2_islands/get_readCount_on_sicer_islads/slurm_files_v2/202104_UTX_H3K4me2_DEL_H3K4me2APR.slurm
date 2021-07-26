#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files_v2/slurm_202104_UTX_H3K4me2_DEL_H3K4me2APR.out

python get_readCount_on_region_from_BamBed.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer2_islands/merged_islands/H3K4me2APR.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/DEL_H3K4me2APR/DEL_H3K4me2APR_treat.bam \
-s hg38 -f bam -o readCount_csv/202104_UTX_H3K4me2/DEL_H3K4me2APR_on_sicer_islands.csv

