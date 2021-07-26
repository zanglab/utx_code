#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_DEL_UTXFEB_on_H3K27ac.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f8_diff_binding_on_ALL_k4me1_k27ac_regions_RPKM/data_hm_peaks/peaks_merged/H3K27ac.merged.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/DEL_UTXFEB/DEL_UTXFEB_treat.bam \
-s hg38 -f bam -m -e 2000 -o rpkm_csv/DEL_UTXFEB_on_H3K27ac_es2kb.csv

