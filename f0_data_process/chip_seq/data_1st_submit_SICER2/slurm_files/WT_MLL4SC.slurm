#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_WT_MLL4SC.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out"

time sicer -t ${processed_bam_dir}/WT_MLL4SC/WT_MLL4SC_treat.bam \
-s hg38 --output_directory sicer_results/WT_MLL4SC
