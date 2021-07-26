#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o run1.out 

module load samtools

for ii in ../processed_out/*/*.sam
do 
 echo $ii
 samtools flagstat -@8 $ii
done