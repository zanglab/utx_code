#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms_k4m1_k27ac/slurm_MT2_H3K27ac.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out//MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -o rpkm_csv/MT2_H3K27ac_on_Promoter.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out//MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -o rpkm_csv/MT2_H3K27ac_on_UDHS.csv
