
data_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq
UTX_peaks=$data_dir/202011_UTX_trim/process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_peaks.narrowPeak
UTX_islands=$data_dir/sicer_df/sicer_out/WT_over_Vector_UTX/WT_UTX_treat-vs-Vector_UTX_treat-W200-G600-E1000-union.island
UTX_islands_increased=$data_dir/sicer_df/sicer_out/WT_over_Vector_UTX/WT_UTX_treat-W200-G600-increased-islands-summary-FDR0.01
UTXFEB_peaks=$data_dir/202102_UTX_H3K27me3_trim/process_qc_out/WT_UTXFEB_with_Vector_control/WT_UTXFEB_with_Vector_control_q01_peaks.narrowPeak
UTXFEB_islands=$data_dir/sicer_df/sicer_out/WT_over_Vector_UTXFEB/WT_UTXFEB_treat-vs-Vector_UTXFEB_treat-W200-G600-E1000-union.island
UTXFEB_islands_increased=$data_dir/sicer_df/sicer_out/WT_over_Vector_UTXFEB/WT_UTXFEB_treat-W200-G600-increased-islands-summary-FDR0.01

cp $UTX_peaks UTX_peaks.bed
cp $UTX_islands UTX_islands.bed
cp $UTX_islands_increased UTX_islands_increased.bed
cp $UTXFEB_peaks UTXFEB_peaks.bed
cp $UTXFEB_islands UTXFEB_islands.bed
cp $UTXFEB_islands_increased UTXFEB_islands_increased.bed

cat UTX_islands.bed UTXFEB_islands.bed > UTX_all_islands.bed
bedtools sort -i UTX_all_islands.bed > UTX_all_islands_sort.bed
bedtools merge -i UTX_all_islands_sort.bed > UTX_all_islands_merge.bed

