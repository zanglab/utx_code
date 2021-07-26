#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_frip/del_cIDR_anti_FLAG.out 

genome_index="/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa"
fqdir="/nv/vol190/zanglab/bs6yk/BS202009_JiangLab/Fastq/"
outdir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202009/process_qc_out/del_cIDR_anti_FLAG"

python 2_chipseq_qc_se_append_add_frip_score.py -i name_list_chipseq_se.txt -o slurm_frip/del_cIDR_anti_FLAG_qc.log -s hg38 -b del_cIDR_anti_FLAG -f .fastq.gz

