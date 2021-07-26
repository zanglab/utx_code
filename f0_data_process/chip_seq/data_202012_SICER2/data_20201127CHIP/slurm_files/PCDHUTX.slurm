#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_PCDHUTX.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_202012/data_20201127CHIP/process_qc_out/PCDHUTX"

time sicer -t ${processed_bam_dir}/PCDHUTX_treat.bam \
-s hg38 --output_directory sicer_results/PCDHUTX
