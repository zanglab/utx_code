#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files_v2/slurm_DEL_over_WT_H3K4me2APR.out

time sicer_df -t \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/DEL_H3K4me2APR/DEL_H3K4me2APR_treat.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/WT_H3K4me2APR/WT_H3K4me2APR_treat.bam \
-s hg38 --output_directory sicer_out/DEL_over_WT_H3K4me2APR
