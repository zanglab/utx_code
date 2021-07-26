#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_pattern_at_PROseq_union_peak_HM_overlapped_slurms/slurm_PCDH2PRO.out

python ../data_modules_revised/get_pattern_near_site_RPKM.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f6_proseq/data_union_peaks/overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/PCDH2PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o f4_pattern_at_PROseq_union_peak_HM_overlapped/PCDH2PRO_PROseq_peak_es2kb_bin200_RPKM.csv

python ../data_modules_revised/get_pattern_near_site_raw_count.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f6_proseq/data_union_peaks/overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/PCDH2PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o f4_pattern_at_PROseq_union_peak_HM_overlapped/PCDH2PRO_PROseq_peak_es2kb_bin200_rawCount.csv

