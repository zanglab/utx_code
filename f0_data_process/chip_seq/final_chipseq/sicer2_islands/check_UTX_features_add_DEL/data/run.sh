cp ../../../202011_UTX_trim/process_qc_out/DEL_UTX_with_Vector_control/DEL_UTX_with_Vector_control_peaks.narrowPeak  .
cp ../../../202011_UTX_trim/process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_peaks.narrowPeak .

cat WT_UTX_with_Vector_control_peaks.narrowPeak |grep -v v > WT_peaks.bed
cat DEL_UTX_with_Vector_control_peaks.narrowPeak |grep -v v > DEL_peaks.bed

bedtools intersect -a DEL_peaks.bed  -b WT_peaks.bed  -wa > DEL_overlap_with_WT.bed
bedtools intersect -a DEL_peaks.bed  -b WT_peaks.bed  -wa -v > DEL_NOT_overlap_with_WT.bed
cat DEL_overlap_with_WT.bed|sort -u > DEL_overlap_with_WT.bed.uniq

bedtools intersect -b DEL_peaks.bed  -a WT_peaks.bed  -wa > WT_overlap_with_DEL.bed
bedtools intersect -b DEL_peaks.bed  -a WT_peaks.bed  -wa -v > WT_NOT_overlap_with_DEL.bed
cat WT_overlap_with_DEL.bed|sort -u > WT_overlap_with_DEL.bed.uniq

# seems bedtools will output repeated lines in intersect	 
