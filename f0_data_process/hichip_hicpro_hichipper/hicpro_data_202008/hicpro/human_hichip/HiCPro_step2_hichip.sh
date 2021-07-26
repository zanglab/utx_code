#!/bin/bash

#SBATCH -n 4
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=30gb
#SBATCH -p parallel
#SBATCH -A cphg_cz3d

#SBATCH --mail-user=zhenjia@virginia.edu
#SBATCH --mail-type=end
#SBATCH --job-name=HiCpro_s1_hichip
#SBATCH --export=ALL

cd $SLURM_SUBMIT_DIR

make --file /nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip_hicpro_hichipper/hicpro/config_mboi.txt CONFIG_SYS=/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/config-system.txt all_persample 2>&1
