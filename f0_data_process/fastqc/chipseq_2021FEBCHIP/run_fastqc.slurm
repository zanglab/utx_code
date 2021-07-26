#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o fastqc.out


file_dir="/nv/vol190/zanglab/bs6yk/2021FEBCHIP"

module load fastqc
mkdir fastqc_out

for fq in ${file_dir}/*L001*.gz
do 
  echo $fq
  fastqc $fq -o fastqc_out 
done




