
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
outdir=overlapped
mkdir $outdir

MLL4_islands=${project_dir}/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_MLL4/WT_MLL4_treat-vs-Vector_MLL4_treat-W200-G600-E1000-union.island

UTX_202011_WT_over_Vector_peaks=${project_dir}/f0_data_process/chip_seq/final_chipseq/202011_UTX_trim/process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_peaks.narrowPeak
UTX_202102_WT_over_Vector_peaks=${project_dir}/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/WT_UTXFEB_with_Vector_control/WT_UTXFEB_with_Vector_control_q01_peaks.narrowPeak

UTX_202011_WT_over_Vector_islands=${project_dir}/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_UTX/WT_UTX_treat-W200-G600-increased-islands-summary-FDR0.01
UTX_202102_WT_over_Vector_islands=${project_dir}/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_UTXFEB/WT_UTXFEB_treat-W200-G600-increased-islands-summary-FDR0.01


# bedtools intersect -a $MLL4_islands -b $UTX_202011_WT_over_Vector_peaks -wa -u > $outdir/MLL4_on_202011UTX_WT_over_Vector_peaks.bed
# bedtools intersect -a $MLL4_islands -b $UTX_202102_WT_over_Vector_peaks -wa -u > $outdir/MLL4_on_202102UTX_WT_over_Vector_peaks.bed
# bedtools intersect -a $MLL4_islands -b $UTX_202011_WT_over_Vector_islands -wa -u > $outdir/MLL4_on_202011UTX_WT_over_Vector_islands.bed
# bedtools intersect -a $MLL4_islands -b $UTX_202102_WT_over_Vector_islands -wa -u > $outdir/MLL4_on_202102UTX_WT_over_Vector_islands.bed

python_file=/nv/vol190/zanglab/zw5j/scripts/Modules/find_overlap_keep_info_NOT_sep_strand.py
python $python_file -a $MLL4_islands -b $UTX_202011_WT_over_Vector_peaks -s hg38 -p $outdir/MLL4_on_202011UTX_WT_over_Vector_peaks.bed
python $python_file -a $MLL4_islands -b $UTX_202102_WT_over_Vector_peaks -s hg38 -p $outdir/MLL4_on_202102UTX_WT_over_Vector_peaks.bed
python $python_file -a $MLL4_islands -b $UTX_202011_WT_over_Vector_islands -s hg38 -p $outdir/MLL4_on_202011UTX_WT_over_Vector_islands.bed
python $python_file -a $MLL4_islands -b $UTX_202102_WT_over_Vector_islands -s hg38 -p $outdir/MLL4_on_202102UTX_WT_over_Vector_islands.bed

