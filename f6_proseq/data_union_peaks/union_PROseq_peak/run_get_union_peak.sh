mkdir peak_files
cp /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_trim/processed_out/*/*narrow* peak_files
cat peak_files/* > proseq_all_peaks.bed
bedtools sort -i proseq_all_peaks.bed > proseq_all_peaks.sorted.bed
bedtools merge -i proseq_all_peaks.sorted.bed > proseq_all_peaks.sorted.merged.bed

cat proseq_all_peaks.sorted.merged.bed |awk '{OFS="\t"; print$0,NR}' >proseq_all_peaks.sorted.merged.id.bed