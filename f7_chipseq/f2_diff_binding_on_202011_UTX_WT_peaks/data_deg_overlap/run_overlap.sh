
peak_file=WT_UTX_with_Vector_control_peaks.narrowPeak

all_genes=gene_promoters/hg38_4k_promoter_geneID.bed
down_genes=gene_promoters/WT_vs_ctrl_Vector_dngenes_promoter.bed
up_genes=gene_promoters/WT_vs_ctrl_Vector_upgenes_promoter.bed

# mkdir peak_overlap_promoter
# bedtools intersect -a $peak_file -b $all_genes -wa -u > peak_overlap_promoter/peak_overlap_ALL_genes.bed
# bedtools intersect -a $peak_file -b $down_genes -wa -u > peak_overlap_promoter/peak_overlap_Up_genes.bed
# bedtools intersect -a $peak_file -b $up_genes -wa -u > peak_overlap_promoter/peak_overlap_Down_genes.bed

for e2 in 0 3 8 98
do
val=`expr $e2 + 2`
echo $val
python find_overlap_keep_info_NOT_sep_strand.py -a $peak_file -b $all_genes -s hg38 -e2 $e2 -p peak_overlap_promoter/peak_overlap_ALL_genes_es${val}kb.bed
python find_overlap_keep_info_NOT_sep_strand.py -a $peak_file -b $down_genes -s hg38 -e2 $e2 -p peak_overlap_promoter/peak_overlap_Down_genes_es${val}kb.bed
python find_overlap_keep_info_NOT_sep_strand.py -a $peak_file -b $up_genes -s hg38 -e2 $e2 -p peak_overlap_promoter/peak_overlap_Up_genes_es${val}kb.bed
done