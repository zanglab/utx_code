#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_at_PROseq_union_peak/slurm_PCDH2PRO.out

python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f6_proseq/data_0x2_MAPQ10/data/union_PROseq_peak/proseq_all_peaks.sorted.merged.id.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/PCDH2PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -o f5_proseq_union_peak_RPKM/PCDH2PRO_on_union_Proseq_peak.csv

python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f6_proseq/data_0x2_MAPQ10/data/union_PROseq_peak/proseq_all_peaks.sorted.merged.id.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10/PCDH2PRO_5prime_strandRev.bed \
-s hg38 -f bed -g 0 -m -e 2000 -o f5_proseq_union_peak_RPKM/PCDH2PRO_on_union_Proseq_peak_es2kb.csv

