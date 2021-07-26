#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_MT2_H3K27ac.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_MLL4/WT_MLL4_treat-vs-Vector_MLL4_treat-W200-G600-E1000-union.island \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -m -e 500 -o rpkm_csv/MT2_H3K27ac_on_MLL4_islands_es500bp.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_MLL4/WT_MLL4_treat-vs-Vector_MLL4_treat-W200-G600-E1000-union.island \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -m -e 2000 -o rpkm_csv/MT2_H3K27ac_on_MLL4_islands_es2kb.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_MLL4/WT_MLL4_treat-vs-Vector_MLL4_treat-W200-G600-E1000-union.island \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -o rpkm_csv/MT2_H3K27ac_on_MLL4_islands.csv

