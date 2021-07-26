#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_EIFUTX4_with_Vector_control.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_4th_20201220/process_qc_out/EIFUTX4_with_Vector_control"

time sicer -t ${processed_bam_dir}/EIFUTX4_with_Vector_control_treat.bam \
-c ${processed_bam_dir}/EIFUTX4_with_Vector_control_control.bam \
-s hg38 --output_directory sicer_results/EIFUTX4_with_Vector_control
