#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o run_split_slurms/slurm_PCDH4HICHIP_1.out

time split_reads.py -n 10000000 --results_folder check_on_terminal_split_fq/PCDH4HICHIP /nv/vol190/zanglab/bs6yk/HICHIP/PCDH4HICHIP_R1.fastq.gz
