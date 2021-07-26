#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o run_split_slurms/slurm_K4M3DEL3_2.out

time split_reads.py -n 20000000 --results_folder check_on_terminal_split_fq/K4M3DEL3 /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/data_20200730/usftp21.novogene.com/raw_data//K4M3DEL3/K4M3DEL3_CKDL200158878-1a_H5T7HDSXY_L3_2.fq.gz
