#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_del_cIDR_UTX_WT_UTX.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_1st_submit/process_qc_out"

time sicer_df -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_1st_submit/process_qc_out/del_cIDR_UTX/del_cIDR_UTX_treat.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_1st_submit/process_qc_out/WT_UTX/WT_UTX_treat.bam \
-s hg38 --output_directory sicer_results/del_cIDR_UTX_vs_WT_UTX
