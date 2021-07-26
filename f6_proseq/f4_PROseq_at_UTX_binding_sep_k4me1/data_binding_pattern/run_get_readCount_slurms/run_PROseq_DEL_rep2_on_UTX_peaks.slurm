#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_readCount_slurms/slurm_PROseq_DEL_rep2_on_UTX_peaks.out

python get_pattern_near_site_readcount_write_out_revised.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f4_strand_specific_bedfile/DEL32PRO_5prime_plus.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o readCount_csv/PROseq_DEL_rep2_plus_on_UTX_peaks_es2kb_bin200.csv

python get_pattern_near_site_readcount_write_out_revised.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f4_strand_specific_bedfile/DEL32PRO_5prime_minus.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o readCount_csv/PROseq_DEL_rep2_minus_on_UTX_peaks_es2kb_bin200.csv

