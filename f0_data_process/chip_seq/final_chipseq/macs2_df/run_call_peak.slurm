#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_macs_peaks.log


mkdir macs2_resluts
data_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/
Vector_H3K27ac=${data_dir}/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K27ac/Vector_H3K27ac_treat.bam
WT_H3K27ac=${data_dir}/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K27ac/WT_H3K27ac_treat.bam

Vector_H3K4me1=${data_dir}/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K4me1/Vector_H3K4me1_treat.bam
WT_H3K4me1=${data_dir}/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K4me1/WT_H3K4me1_treat.bam

Vector_MLL4=${data_dir}/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out/Vector_MLL4/Vector_MLL4_treat.bam
WT_MLL4=${data_dir}/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out/WT_MLL4/WT_MLL4_treat.bam


## ==== H3K27ac
macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
$WT_H3K27ac -c \
$Vector_H3K27ac -n macs2_resluts/H3K27ac_WT_with_Vector_control

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
$WT_H3K27ac -c \
$Vector_H3K27ac -n macs2_resluts/H3K27ac_WT_with_Vector_control_q01


## ==== H3K4me1
macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
$WT_H3K4me1 -c \
$Vector_H3K4me1 -n macs2_resluts/H3K4me1_WT_with_Vector_control

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
$WT_H3K4me1 -c \
$Vector_H3K4me1 -n macs2_resluts/H3K4me1_WT_with_Vector_control_q01


## ==== MLL4
macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
$WT_MLL4 -c \
$Vector_MLL4 -n macs2_resluts/MLL4_WT_with_Vector_control

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
$WT_MLL4 -c \
$Vector_MLL4 -n macs2_resluts/MLL4_WT_with_Vector_control_q01



