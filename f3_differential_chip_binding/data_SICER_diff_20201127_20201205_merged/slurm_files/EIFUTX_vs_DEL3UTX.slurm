#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_EIFUTX_DEL3UTX.out

processed_bam_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_202012/data_20201127CHIP_20201205CHIP_merged/process_qc_out"

time sicer_df -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_202012/data_20201127CHIP_20201205CHIP_merged/process_qc_out/EIFUTX/EIFUTX_treat.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq//data_202012/data_20201127CHIP_20201205CHIP_merged/process_qc_out/DEL3UTX/DEL3UTX_treat.bam \
-s hg38 --output_directory sicer_results/EIFUTX_vs_DEL3UTX
