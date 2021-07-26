#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 1:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_py4.log

# Run program

time python py4_heatmap_across_CellTypes.py
