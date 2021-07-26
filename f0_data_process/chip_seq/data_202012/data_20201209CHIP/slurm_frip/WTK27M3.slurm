#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_frip/WTK27M3.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/Jianglab202012CHIP/20201209CHIP/"
outdir="./process_qc_out/WTK27M3"

python 2_chipseq_qc_se_append_add_frip_score.py -i name_list_chipseq_se.txt -o slurm_frip/WTK27M3_qc.log -s hg38 -b WTK27M3 -f .fastq.gz

