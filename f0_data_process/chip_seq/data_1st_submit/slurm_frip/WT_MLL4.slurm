#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_frip/WT_MLL4.out 

python 2_chipseq_qc_se_append_add_frip_score.py -i name_list_chipseq_se_append.txt -o slurm_frip/WT_MLL4_qc.log -s hg38 -b WT_MLL4 -f .fastq.gz

