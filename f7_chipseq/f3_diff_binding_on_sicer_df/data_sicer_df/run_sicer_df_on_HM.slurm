#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_sider_df.out


vector_k4=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K4me1/Vector_H3K4me1_treat.bam
vector_k27=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K27ac/Vector_H3K27ac_treat.bam
wt_k4=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K4me1/WT_H3K4me1_treat.bam
wt_k27=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K27ac/WT_H3K27ac_treat.bam

outdir=sicer_df_results
mkdir $outdir

time sicer_df -t $wt_k4 $vector_k4 -s hg38 \
--output_directory ${outdir}/WT_vs_Vector_H3K4me1

time sicer_df -t $wt_k27 $vector_k27 -s hg38 \
--output_directory ${outdir}/WT_vs_Vector_H3K27ac
