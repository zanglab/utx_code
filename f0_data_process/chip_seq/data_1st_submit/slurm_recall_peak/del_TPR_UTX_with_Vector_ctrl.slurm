#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_recall_peak/del_TPR_UTX_with_Vector_ctrl.out 

macs2 callpeak --SPMR -B -q 0.2 --keep-dup 1 -g hs -t \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_DTPR-UTX_L001.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_DTPR-UTX_L002.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_DTPR-UTX_L003.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_DTPR-UTX_L004.bam -c \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_PCDH-UTX_L001.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_PCDH-UTX_L002.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_PCDH-UTX_L003.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/del_TPR_UTX_with_Vector_ctrl_PCDH-UTX_L004.bam -n /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/del_TPR_UTX_with_Vector_ctrl/recall_peak_q02_del_TPR_UTX_with_Vector_ctrl

