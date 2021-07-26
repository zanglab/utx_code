#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_macs_peaks_UTX.log


#mkdir macs2_resluts
data_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/

Vector_UTX=${data_dir}/202011_UTX_trim/process_qc_out/Vector_UTX/Vector_UTX_treat.bam
WT_UTX=${data_dir}/202011_UTX_trim/process_qc_out/WT_UTX/WT_UTX_treat.bam
DEL_UTX=${data_dir}/202011_UTX_trim/process_qc_out/DEL_UTX/DEL_UTX_treat.bam

## ==== UTX
macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
$WT_UTX -c \
$Vector_UTX -n macs2_resluts/UTX_WT_with_Vector_control

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
$DEL_UTX -c \
$Vector_UTX -n macs2_resluts/UTX_DEL_with_Vector_control


