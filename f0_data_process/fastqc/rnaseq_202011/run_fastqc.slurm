#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o fastqc.out


file_dir="/nv/vol190/zanglab/bs6yk/202011RNASEQ//"

mkdir fastqc_out
module load fastqc

for fq in ${file_dir}/*L001*.gz
do 
  echo $fq
  fastqc $fq -o fastqc_out 
done




