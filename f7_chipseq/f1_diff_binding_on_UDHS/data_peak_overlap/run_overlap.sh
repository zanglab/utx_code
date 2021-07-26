
promoter_file='/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'
udhs_file='/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
mkdir peak_files
mkdir merged_peaks
mkdir overlapped

## ==== copy the peak files
peak_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang//f0_data_process/chip_seq/final_chipseq
cp ${peak_dir}/*trim/process_qc_out/*/*narrow*  peak_files
rm peak_files/*_q01_peaks.narrowPeak


# ==== merge all peaks of each factor and check overlap with UDHS
for antibody in UTX MLL4 H3K4me1 H3K4me2 H3K4me3 H3K27ac
do
  #echo peak_files/*_${antibody}_*
  cat peak_files/*_${antibody}_*  > merged_peaks/${antibody}_all_peaks.bed
  bedtools sort -i merged_peaks/${antibody}_all_peaks.bed > merged_peaks/${antibody}_all_peaks.sorted.bed
  bedtools merge -i merged_peaks/${antibody}_all_peaks.sorted.bed > merged_peaks/${antibody}_all_peaks.merged.bed
  bedtools intersect -a $promoter_file -b merged_peaks/${antibody}_all_peaks.merged.bed -wa -u > overlapped/Promoter_overlap_${antibody}.bed
  bedtools intersect -a $udhs_file -b merged_peaks/${antibody}_all_peaks.merged.bed -wa -u > overlapped/UDHS_overlap_${antibody}.bed
done


# ==== merge all peaks of ALL factorand check overlap with UDHS
cat peak_files/*narrowPeak  > merged_peaks/ALL_peaks.bed
bedtools sort -i merged_peaks/ALL_peaks.bed > merged_peaks/ALL_peaks.sorted.bed
bedtools merge -i merged_peaks/ALL_peaks.sorted.bed > merged_peaks/ALL_peaks.merged.bed
bedtools intersect -a $promoter_file -b merged_peaks/ALL_peaks.merged.bed -wa -u > overlapped/Promoter_overlap_ALL.bed
bedtools intersect -a $udhs_file -b merged_peaks/ALL_peaks.merged.bed -wa -u > overlapped/UDHS_overlap_ALL.bed


## ==== merged peak on UDHS
mkdir peaks_on_UDHS
for ii in merged_peaks/*.merged.bed
do
  outname=$(basename $ii .bed)
  bedtools intersect  -a $ii -b $udhs_file -wa -u > peaks_on_UDHS/${outname}_on_UDHS.bed
done


