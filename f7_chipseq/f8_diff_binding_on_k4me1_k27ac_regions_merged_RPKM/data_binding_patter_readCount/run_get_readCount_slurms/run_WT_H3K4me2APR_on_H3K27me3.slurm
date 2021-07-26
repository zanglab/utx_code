#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_readCount_slurms/slurm_WT_H3K4me2APR_on_H3K27me3.out

python get_pattern_near_site_readcount_write_out_revised.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f8_diff_binding_on_ALL_k4me1_k27ac_regions_RPKM/data_hm_peaks/peaks_merged/H3K27me3.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/WT_H3K4me2APR/WT_H3K4me2APR_treat.bam \
-s hg38 -f bam -w 2000 -b 200 -m -o readCount_csv/WT_H3K4me2APR_on_H3K27me3_es2kb_bin200.csv

