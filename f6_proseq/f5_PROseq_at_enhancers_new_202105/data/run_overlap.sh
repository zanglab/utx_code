
bedtools intersect \
-a union_PROseq_peak/proseq_all_peaks.sorted.merged.id.bed \
-b overlapped/H3K4me1_overlapping_H3K27ac.islands.bed \
-wa -u > overlapped/union_PROseq_peak_overlapping_H3K4me1_H3K27ac.bed


mkdir overlapped

### PROseq overlapping H3K4me1_overlapping_H3K27ac w/ UTX
python find_overlap_keep_info_NOT_sep_strand.py \
-a H3K4me1_overlapping_H3K27ac.islands.bed \
-b WT_UTX_with_Vector_control_peaks.narrowPeak \
-s hg38 -p overlapped/H3K4me1_overlapping_H3K27ac_overlapping_UTX.bed

python find_overlap_keep_info_NOT_sep_strand.py \
-a proseq_all_peaks.sorted.merged.id.bed \
-b overlapped/H3K4me1_overlapping_H3K27ac_overlapping_UTX.bed \
-s hg38 -p overlapped/PROseq_on_H3K4me1_overlapping_H3K27ac_overlapping_UTX.bed


### PROseq overlapping UTX that are on H3K4me1_overlapping_H3K27ac
python find_overlap_keep_info_NOT_sep_strand.py \
-a WT_UTX_with_Vector_control_peaks.narrowPeak \
-b H3K4me1_overlapping_H3K27ac.islands.bed \
-s hg38 -p overlapped/UTX_on_H3K4me1_overlapping_H3K27ac.bed

python find_overlap_keep_info_NOT_sep_strand.py \
-a proseq_all_peaks.sorted.merged.id.bed \
-b overlapped/UTX_on_H3K4me1_overlapping_H3K27ac.bed \
-s hg38 -p overlapped/PROseq_on_UTX_on_H3K4me1_overlapping_H3K27ac.bed



## PROseq overlapping either H3K27ac w/ UTX or H3K4me1 w/ UTX
python find_overlap_keep_info_NOT_sep_strand.py \
-a H3K27ac.merged.bed \
-b WT_UTX_with_Vector_control_peaks.narrowPeak \
-s hg38 -p overlapped/H3K27ac_with_UTX.bed

python find_overlap_keep_info_NOT_sep_strand.py \
-a H3K4me1.merged.bed \
-b WT_UTX_with_Vector_control_peaks.narrowPeak \
-s hg38 -p overlapped/H3K4me1_with_UTX.bed

cat overlapped/H3K27ac_with_UTX.bed overlapped/H3K4me1_with_UTX.bed > overlapped/H3K27ac_H3K4me1_with_UTX.bed

python find_overlap_keep_info_NOT_sep_strand.py \
-a proseq_all_peaks.sorted.merged.id.bed \
-b overlapped/H3K27ac_H3K4me1_with_UTX.bed \
-s hg38 -p overlapped/PROseq_on_H3K27ac_H3K4me1_with_UTX.bed


## for those PRO-seq peaks at enchancers
## separate into promoter overlapped and promoter NON overlapped groups
promoter_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed.tss.sort.merge.bed
python find_overlap_keep_info_NOT_sep_strand.py \
-a overlapped/PROseq_on_H3K27ac_H3K4me1_with_UTX.bed \
-b $promoter_file \
-s hg38 \
-p overlapped/PROseq_on_H3K27ac_H3K4me1_with_UTX_promoter_overlapped.bed \
-q overlapped/PROseq_on_H3K27ac_H3K4me1_with_UTX_promoter_NON_overlapped.bed









