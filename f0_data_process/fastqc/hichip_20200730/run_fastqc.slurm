#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o fastqc.out


file_dir="/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/data_20200730/usftp21.novogene.com/raw_data"

mkdir fastqc_out

for fq in ${file_dir}/*/*1.fq.gz
do 
  echo $fq
  fastqc $fq -o fastqc_out 
done




