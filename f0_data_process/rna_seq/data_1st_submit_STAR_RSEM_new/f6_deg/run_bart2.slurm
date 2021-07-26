#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o out_bart2.log



for ii in f2_deg/*.txt
do 
    bart2 geneset -i $ii -s hg38 --outdir f3_deg_bart2
done
