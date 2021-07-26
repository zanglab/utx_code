#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_UTX_eIFIDR_MLL4SC.out

python get_RPM_RPKM_count_on_regions.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_MLL4SC/UTX_eIFIDR_MLL4SC_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_MLL4SC_on_Promoter.csv
python get_RPM_RPKM_count_on_regions.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_MLL4SC/UTX_eIFIDR_MLL4SC_treat.bam -s hg38 -f bam -m -e 10000 -o rpkm_csv/UTX_eIFIDR_MLL4SC_on_Promoter_es10kb.csv
python get_RPM_RPKM_count_on_regions.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_unionDHS_fc4_50merge.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_MLL4SC/UTX_eIFIDR_MLL4SC_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_MLL4SC_on_UDHS.csv
