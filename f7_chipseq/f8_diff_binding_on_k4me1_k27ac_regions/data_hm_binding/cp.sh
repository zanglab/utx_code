
data_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq

H3K27ac_islands_increased=$data_dir/sicer_df/sicer_out/WT_over_Vector_H3K27ac/WT_H3K27ac_treat-W200-G600-increased-islands-summary-FDR0.01
H3K4me1_islands_increased=$data_dir/sicer_df/sicer_out/WT_over_Vector_H3K4me1/WT_H3K4me1_treat-W200-G600-increased-islands-summary-FDR0.01
MLL4_islands_increased=$data_dir/sicer_df/sicer_out/WT_over_Vector_MLL4/WT_MLL4_treat-W200-G600-increased-islands-summary-FDR0.01

cp $H3K27ac_islands_increased H3K27ac_islands_increased.bed
cp $H3K4me1_islands_increased H3K4me1_islands_increased.bed
cp $MLL4_islands_increased MLL4_islands_increased.bed


H3K27ac_peaks=$data_dir/macs2_df/macs2_resluts/H3K27ac_WT_with_Vector_control_q01_peaks.narrowPeak
H3K4me1_peaks=$data_dir/macs2_df/macs2_resluts/H3K4me1_WT_with_Vector_control_q01_peaks.narrowPeak
MLL4_peaks=$data_dir/macs2_df/macs2_resluts/MLL4_WT_with_Vector_control_q01_peaks.narrowPeak

cp $H3K27ac_peaks H3K27ac_peaks.bed
cp $H3K4me1_peaks H3K4me1_peaks.bed
cp $MLL4_peaks MLL4_peaks.bed

