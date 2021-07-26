#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_readCount_slurms/slurm_WT_UTXFEB_on_UTXFEB_peaks.out

python get_pattern_near_site_readcount_write_out_revised.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTXFEB_peaks.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/WT_UTXFEB/WT_UTXFEB_treat.bam \
-s hg38 -f bam -w 2000 -b 200 -m -o readCount_csv/WT_UTXFEB_on_UTXFEB_peaks_es2kb_bin200.csv

