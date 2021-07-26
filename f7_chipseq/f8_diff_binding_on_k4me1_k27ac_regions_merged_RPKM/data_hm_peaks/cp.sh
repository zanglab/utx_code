
outdir=peaks_cp
mkdir $outdir

#copy the hm peaks
data_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq
H3K4me1_peaks=$data_dir/202102_H3K27ac_H3K4me1_trim/process_qc_out/*H3K4me1/*.narrowPeak
H3K27ac_peaks=$data_dir/202102_H3K27ac_H3K4me1_trim/process_qc_out/*H3K27ac/*.narrowPeak
H3K4me3_peaks=$data_dir/re_1st_submission_H3K4me3_MLL4SC_trim/process_qc_out/*H3K4me3/*.narrowPeak
H3K27me3_peaks=$data_dir/202102_UTX_H3K27me3_trim/process_qc_out/*H3K27me3/*.narrowPeak

cp $H3K4me1_peaks $outdir
cp $H3K27ac_peaks $outdir
cp $H3K4me3_peaks $outdir
cp $H3K27me3_peaks $outdir

rm $outdir/*_q01_*
rm $outdir/EIF*
rm $outdir/TPR*

#then merge the hm peaks
merge_dir=peaks_merged
mkdir $merge_dir

for hm in H3K4me1 H3K4me3 H3K27ac H3K27me3
do
  cat $outdir/*${hm}*.narrowPeak > $merge_dir/${hm}.all.bed
  bedtools sort -i $merge_dir/${hm}.all.bed > $merge_dir/${hm}.sort.bed
  bedtools merge -i $merge_dir/${hm}.sort.bed > $merge_dir/${hm}.merged.bed
done


