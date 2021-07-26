# copy the union of peaks from each factor
mkdir union_factor_peaks 
cp /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f1_differential_binding_on_UDHS/data_peak_overlap/merged_peaks/UTX*.merged.bed union_factor_peaks/
cp /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/sicer2_islands/merged_islands/H*.merge.bed union_factor_peaks/
cp /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f1_differential_binding_on_UDHS/data_peak_overlap/peak_files/WT_UTX_with_Vector_control_peaks.narrowPeak union_factor_peaks/

# potential enhancers with K4, K27 and PROseq peaks
mkdir overlapped
bedtools intersect \
-a union_factor_peaks/H3K4me1.merge.bed \
-b union_factor_peaks/H3K27ac.merge.bed \
-wa -u > overlapped/H3K4me1_overlapping_H3K27ac.islands.bed

bedtools intersect \
-a union_PROseq_peak/proseq_all_peaks.sorted.merged.id.bed \
-b overlapped/H3K4me1_overlapping_H3K27ac.islands.bed \
-wa -u > overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed

# separate potential enhancers into UTX-bound and UTX-unbound regions
bedtools intersect \
-a overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed \
-b union_factor_peaks/WT_UTX_with_Vector_control_peaks.narrowPeak \
-wa -u > overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac_with_WT_UTX_bound.bed

bedtools intersect \
-a overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed \
-b union_factor_peaks/UTX_all_peaks.merged.bed \
-wa -v > overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac_with_ALL_UTX_unbound.bed
