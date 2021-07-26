#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm.out



time python py1_get_peak_overlapped_UDHS.py 