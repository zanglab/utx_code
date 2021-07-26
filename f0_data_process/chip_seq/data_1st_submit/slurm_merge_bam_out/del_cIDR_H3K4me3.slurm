#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_merge_bam_out/del_cIDR_H3K4me3.out 

samtools merge -@ 4 /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_H3K4me3/del_cIDR_H3K4me3_treat.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_H3K4me3/del_cIDR_H3K4me3_Del3-K4M3_L001.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_H3K4me3/del_cIDR_H3K4me3_Del3-K4M3_L002.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_H3K4me3/del_cIDR_H3K4me3_Del3-K4M3_L003.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_cIDR_H3K4me3/del_cIDR_H3K4me3_Del3-K4M3_L004.bam

