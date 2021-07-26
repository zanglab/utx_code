#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_readCount_slurms/slurm_MT2_H3K27ac.out

python get_pattern_near_site_readcount_write_out.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang//f0_data_process/chip_seq/final_chipseq/202011_UTX_trim/process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_peaks.narrowPeak \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/MT2_H3K27ac/MT2_H3K27ac_treat.bam \
-s hg38 -f bam -w 2000 -b 200 -m -o readCount_csv/MT2_H3K27ac_es2kb_bin200_on_202011_UTX_WT_peaks.csv

