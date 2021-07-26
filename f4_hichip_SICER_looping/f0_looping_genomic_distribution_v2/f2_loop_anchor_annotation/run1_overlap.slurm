#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run1.log
#Run program


python_script="find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py"
promoter_file="data_genomic/hg38_4k_promoter_geneID.bed.sorted.merged"
exon_file="data_genomic/hg38_exons.bed.sorted.merged"
intron_file="data_genomic/hg38_introns.bed.sorted.merged"
gb_file="data_genomic/hg38_minus2kb_and_GB_geneID.bed.sorted.merged"
udhs_file="data_genomic/hg38_unionDHS_fc4_50merge.bed"


#for subdir in data_1st_submission_rep_combined  data_1st_submission_sep_rep  data_202008
for subdir in data_202008
do
  outdir=f1_overlapped/$subdir
  mkdir -p $outdir
  for bedfile in f0_data_loop_bedpe/$subdir/*
  do
    prename=$(basename $bedfile .bed)
    echo $outdir $prename
    time python $python_script -a $bedfile -b $promoter_file -s hg38 -p $outdir/${prename}_promoter.bed
    time python $python_script -a $bedfile -b $udhs_file     -s hg38 -p $outdir/${prename}_UDHS.bed
    #time python $python_script -a $bedfile -b $exon_file     -s hg38 -p $outdir/${prename}_exon.bed
    #time python $python_script -a $bedfile -b $intron_file   -s hg38 -p $outdir/${prename}_intron.bed
    #time python $python_script -a $bedfile -b $gb_file       -s hg38 -p $outdir/${prename}_GB.bed
  done
done



