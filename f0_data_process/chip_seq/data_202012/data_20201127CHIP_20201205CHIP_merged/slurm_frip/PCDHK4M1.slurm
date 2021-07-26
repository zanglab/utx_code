#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_frip/PCDHK4M1.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/chipseq/chipseq202012/"
outdir="./process_qc_out/PCDHK4M1"

python 2_chipseq_qc_se_append_add_frip_score.py -i name_list_chipseq_se.txt -o slurm_frip/PCDHK4M1_qc.log -s hg38 -b PCDHK4M1 -f .fastq.gz

