#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=30000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_pattern_slurms/slurm_del_cIDR_UTX.out

python get_pattern_near_site_readcount_write_out.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_UTX/del_cIDR_UTX_treat.bam -s hg38 -f bam -w 10000 -b 1000 -m -o pattern_csv/del_cIDR_UTX_TSS_es10kb_bin1000.csv
