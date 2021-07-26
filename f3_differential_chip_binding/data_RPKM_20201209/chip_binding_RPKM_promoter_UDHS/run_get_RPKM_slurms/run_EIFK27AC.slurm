#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_EIFK27AC.out

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202012/data_20201209CHIP/process_qc_out//EIFK27AC/EIFK27AC_treat.bam \
-s hg38 -f bam -o rpkm_csv/EIFK27AC_on_Promoter.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202012/data_20201209CHIP/process_qc_out//EIFK27AC/EIFK27AC_treat.bam \
-s hg38 -f bam -m -e 10000 -o rpkm_csv/EIFK27AC_on_Promoter_es10kb.csv

python get_RPM_RPKM_count_on_regions.py \
-i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_unionDHS_fc4_50merge.bed \
-t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_202012/data_20201209CHIP/process_qc_out//EIFK27AC/EIFK27AC_treat.bam \
-s hg38 -f bam -o rpkm_csv/EIFK27AC_on_UDHS.csv
