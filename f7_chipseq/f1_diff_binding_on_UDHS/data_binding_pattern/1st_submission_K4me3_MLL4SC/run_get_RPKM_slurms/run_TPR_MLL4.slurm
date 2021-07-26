#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_TPR_MLL4.out

python get_pattern_near_site_readcount_write_out.py \
-i /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out//TPR_MLL4/TPR_MLL4_treat.bam \
-s hg38 -f bam -w 2000 -b 200 -m -o rpkm_csv/TPR_MLL4_es2kb_bin200_on_Promoter.csv

python get_pattern_near_site_readcount_write_out.py \
-i /nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out//TPR_MLL4/TPR_MLL4_treat.bam \
-s hg38 -f bam -w 2000 -b 200 -m -o rpkm_csv/TPR_MLL4_es2kb_bin200_on_UDHS.csv

