#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_get_RPKM_slurms/slurm_UTX_eIFIDR_H3K4me1_rep.out

python get_RPKM_and_count_on_regions_readcount.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/macs2_peaks_merged/union_H3K4me1_rep_peaks.narrowPeak.sorted.merged.id -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_H3K4me1_rep_on_Macs2peak.csv
python get_RPKM_and_count_on_regions_readcount.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/sicer2_islands_merged/union_H3K4me1_rep_treat.scoreisland.sorted.merged.id -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_H3K4me1_rep_on_Sicer2peak.csv
python get_RPKM_and_count_on_regions_readcount.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_H3K4me1_rep_on_Promoter.csv
python get_RPKM_and_count_on_regions_readcount.py -i /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_minus2kb_and_GB_geneID.bed -t /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_treat.bam -s hg38 -f bam -o rpkm_csv/UTX_eIFIDR_H3K4me1_rep_on_GeneBody.csv
