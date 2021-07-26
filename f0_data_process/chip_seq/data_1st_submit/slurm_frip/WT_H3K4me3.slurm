#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_frip/WT_H3K4me3.out 

python 2_chipseq_qc_se_append_add_frip_score.py -i name_list_chipseq_se.txt -o slurm_frip/WT_H3K4me3_qc.log -s hg38 -b WT_H3K4me3 -f .fastq.gz

