#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run1_overlap.log
#Run program


# anchor1: left anchor, anchor2: right anchor
# mkdir overlapped
# mkdir overlapped_bedtools_test

python_script="find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py"
promoter_file="data_genomic/hg38_4k_promoter_geneID.bed.sorted.merged"
exon_file="data_genomic/hg38_exons.bed.sorted.merged"
intron_file="data_genomic/hg38_introns.bed.sorted.merged"
gb_file="data_genomic/hg38_minus2kb_and_GB_geneID.bed.sorted.merged"
udhs_file="data_genomic/hg38_unionDHS_fc4_50merge.bed"

for hm in K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT
# for hm in K27ACDEL
do
#     time bedtools intersect -a data_loop_bedpe/${hm}_anchor1.bed -b $promoter_file -wb -loj > overlapped_bedtools_test/${hm}_promoter_anchor1.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor1.bed -b $promoter_file -s hg38 -p overlapped/${hm}_promoter_anchor1.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor1.bed -b $exon_file -s hg38 -p overlapped/${hm}_exon_anchor1.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor1.bed -b $intron_file -s hg38 -p overlapped/${hm}_intron_anchor1.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor1.bed -b $gb_file -s hg38 -p overlapped/${hm}_genebody_anchor1.bed
#     
#     time python $python_script -a data_loop_bedpe/${hm}_anchor2.bed -b $promoter_file -s hg38 -p overlapped/${hm}_promoter_anchor2.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor2.bed -b $exon_file -s hg38 -p overlapped/${hm}_exon_anchor2.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor2.bed -b $intron_file -s hg38 -p overlapped/${hm}_intron_anchor2.bed
#     time python $python_script -a data_loop_bedpe/${hm}_anchor2.bed -b $gb_file -s hg38 -p overlapped/${hm}_genebody_anchor2.bed

    time python $python_script -a data_loop_bedpe/${hm}_anchor1.bed -b $udhs_file -s hg38 -p overlapped/${hm}_udhs_anchor1.bed
    time python $python_script -a data_loop_bedpe/${hm}_anchor2.bed -b $udhs_file -s hg38 -p overlapped/${hm}_udhs_anchor2.bed


done

