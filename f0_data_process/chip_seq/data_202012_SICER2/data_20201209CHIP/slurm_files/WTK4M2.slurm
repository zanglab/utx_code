#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_WTK4M2.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_202012/data_20201209CHIP/process_qc_out/WTK4M2"

time sicer -t ${processed_bam_dir}/WTK4M2_treat.bam \
-s hg38 --output_directory sicer_results/WTK4M2
