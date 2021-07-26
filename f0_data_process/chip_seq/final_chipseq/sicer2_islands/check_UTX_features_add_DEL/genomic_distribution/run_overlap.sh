
mkdir overlapped

data_dir=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/

bedtools intersect -a ../data/DEL_NOT_overlap_with_WT.bed -b ${data_dir}/hg38_exons.bed -wa > overlapped/DEL_specific_exons_overlapped.bed
bedtools intersect -a ../data/DEL_NOT_overlap_with_WT.bed -b ${data_dir}/hg38_introns.bed -wa > overlapped/DEL_specific_introns_overlapped.bed
bedtools intersect -a ../data/DEL_NOT_overlap_with_WT.bed -b ${data_dir}/hg38_4k_promoter_geneID.bed -wa > overlapped/DEL_specific_promoter_overlapped.bed


bedtools intersect -a ../data/WT_peaks.bed -b ${data_dir}/hg38_exons.bed -wa > overlapped/WT_exons_overlapped.bed
bedtools intersect -a ../data/WT_peaks.bed -b ${data_dir}/hg38_introns.bed -wa > overlapped/WT_introns_overlapped.bed
bedtools intersect -a ../data/WT_peaks.bed -b ${data_dir}/hg38_4k_promoter_geneID.bed -wa > overlapped/WT_promoter_overlapped.bed