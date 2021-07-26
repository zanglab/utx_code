#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_EIF_H3K4me2.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_202012_H3K4me2_trim/process_qc_out//EIF_H3K4me2/EIF_H3K4me2_treat.bam \
-s hg38 -f bam -o rpkm_csv/EIF_H3K4me2_on_Promoter.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_202012_H3K4me2_trim/process_qc_out//EIF_H3K4me2/EIF_H3K4me2_treat.bam \
-s hg38 -f bam -o rpkm_csv/EIF_H3K4me2_on_UDHS.csv
