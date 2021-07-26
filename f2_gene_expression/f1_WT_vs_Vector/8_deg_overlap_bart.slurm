#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 48:00:00
#SBATCH -A cphg_cz3d
#SBATCH -p standard
#SBATCH -o out_run8.log 



for genelist in f7_deg_venn/*.txt
do
    bart2 geneset -i $genelist -s hg38 --outdir f8_bart2_results
done


