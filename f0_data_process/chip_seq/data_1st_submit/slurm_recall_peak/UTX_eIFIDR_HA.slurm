#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_recall_peak/UTX_eIFIDR_HA.out 

macs2 callpeak --SPMR -B -q 0.2 --keep-dup 1 -g hs -t \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_HA/UTX_eIFIDR_HA_EIF-HA_L001.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_HA/UTX_eIFIDR_HA_EIF-HA_L002.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_HA/UTX_eIFIDR_HA_EIF-HA_L003.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_HA/UTX_eIFIDR_HA_EIF-HA_L004.bam -n /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_HA/recall_peak_q02_UTX_eIFIDR_HA

