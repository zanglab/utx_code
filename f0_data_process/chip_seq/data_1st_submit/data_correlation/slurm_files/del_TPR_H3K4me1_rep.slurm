#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_del_TPR_H3K4me1_rep.out

time python get_RPKM_on_regions_readcount.py -s hg38 -f bam -i hg38_5000_ord.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_H3K4me1_rep/del_TPR_H3K4me1_rep_treat.bam \
-o RPKM_csv/del_TPR_H3K4me1_rep_RPKM_genome_5kb.csv

time python get_RPKM_on_regions_readcount.py -s hg38 -f bam -i hg38_4k_promoter_refGene.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_H3K4me1_rep/del_TPR_H3K4me1_rep_treat.bam \
-o RPKM_csv/del_TPR_H3K4me1_rep_RPKM_promoter.csv