#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_Vector_H3K27ac_on_H3K27me3.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f8_diff_binding_on_ALL_k4me1_k27ac_regions_RPKM/data_hm_peaks/peaks_merged/H3K27me3.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K27ac/Vector_H3K27ac_treat.bam \
-s hg38 -f bam -m -e 2000 -o rpkm_csv/Vector_H3K27ac_on_H3K27me3_es2kb.csv

