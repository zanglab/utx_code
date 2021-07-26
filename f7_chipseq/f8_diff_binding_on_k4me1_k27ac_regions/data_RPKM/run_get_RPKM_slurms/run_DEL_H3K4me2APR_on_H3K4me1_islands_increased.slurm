#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_DEL_H3K4me2APR_on_H3K4me1_islands_increased.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f8_diff_binding_on_k4me1_k27ac_regions/data_hm_binding/H3K4me1_islands_increased.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/DEL_H3K4me2APR/DEL_H3K4me2APR_treat.bam \
-s hg38 -f bam -m -e 500 -o rpkm_csv/DEL_H3K4me2APR_on_H3K4me1_islands_increased_es500bp.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f8_diff_binding_on_k4me1_k27ac_regions/data_hm_binding/H3K4me1_islands_increased.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/DEL_H3K4me2APR/DEL_H3K4me2APR_treat.bam \
-s hg38 -f bam -m -e 2000 -o rpkm_csv/DEL_H3K4me2APR_on_H3K4me1_islands_increased_es2kb.csv

