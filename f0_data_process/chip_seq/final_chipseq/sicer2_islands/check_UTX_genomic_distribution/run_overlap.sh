
data_dir=/Volumes/zanglab/zw5j/data/geneID_annotation/hg38/

python find_overlap_keep_info_NOT_sep_strand.py -a WT_UTX_with_Vector_control_peaks.narrowPeak -b ${data_dir}/hg38_exons.bed -s hg38 -p UTX_exons_overlapped.bed
python find_overlap_keep_info_NOT_sep_strand.py -a WT_UTX_with_Vector_control_peaks.narrowPeak -b ${data_dir}/hg38_introns.bed -s hg38 -p UTX_introns_overlapped.bed
python find_overlap_keep_info_NOT_sep_strand.py -a WT_UTX_with_Vector_control_peaks.narrowPeak -b ${data_dir}/hg38_4k_promoter_geneID.bed -s hg38 -p UTX_promoter_overlapped.bed

