#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_promoter_GB_UDHS_count_slurms/slurm_FL1PRO.out

python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data//geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f2_processed_data_with_dup_sep_rep/FL1PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -o f3_promoter_GB_UDHS_count/FL1PRO_on_promoter.csv

python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data//geneID_annotation/hg38/hg38_plus500bp_and_GB_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f2_processed_data_with_dup_sep_rep/FL1PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -o f3_promoter_GB_UDHS_count/FL1PRO_on_GB.csv

python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data//unionDHS/hg38_unionDHS_fc4_50merge.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f2_processed_data_with_dup_sep_rep/FL1PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -o f3_promoter_GB_UDHS_count/FL1PRO_on_UDHS.csv

