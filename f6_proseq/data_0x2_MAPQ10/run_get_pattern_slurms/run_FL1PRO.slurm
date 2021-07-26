#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_pattern_slurms/slurm_FL1PRO.out

python ../data_modules_revised/get_pattern_near_site_RPKM.py \
-i /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/FL1PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o f1_TSS_pattern/FL1PRO_TSS_es2kb_bin200_RPKM.csv

python ../data_modules_revised/get_pattern_near_site_raw_count.py \
-i /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/FL1PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o f1_TSS_pattern/FL1PRO_TSS_es2kb_bin200_rawCount.csv

