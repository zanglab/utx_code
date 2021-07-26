#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_WT_UTX_on_UTXFEB_islands.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTXFEB_islands.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202011_UTX_trim/process_qc_out/WT_UTX/WT_UTX_treat.bam \
-s hg38 -f bam -m -e 500 -o rpkm_csv/WT_UTX_on_UTXFEB_islands_es500bp.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTXFEB_islands.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202011_UTX_trim/process_qc_out/WT_UTX/WT_UTX_treat.bam \
-s hg38 -f bam -m -e 2000 -o rpkm_csv/WT_UTX_on_UTXFEB_islands_es2kb.csv

