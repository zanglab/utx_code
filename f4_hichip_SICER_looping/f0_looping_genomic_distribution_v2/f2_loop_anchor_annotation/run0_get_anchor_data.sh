
# ===========
# batch1 data
# ===========

outdir=f0_data_loop_bedpe/data_1st_submission_sep_rep
mkdir -p $outdir
maps_results_dir=../../f0_looping_data_v2/data_1st_submission_reindex

# anchor1: left anchor, anchor2: right anchor
for bedpe_file in $maps_results_dir/*bedpe
do
    outname=$(basename $bedpe_file .bedpe)
    echo $outname
    tail -n +2 $bedpe_file | awk '{OFS="\t";print$2,$3,$4,$1}' > $outdir/${outname}_anchor1.bed
    tail -n +2 $bedpe_file | awk '{OFS="\t";print$5,$6,$7,$1}' > $outdir/${outname}_anchor2.bed
done

# ==== then merge the data in two replicates
indir=f0_data_loop_bedpe/data_1st_submission_sep_rep
outdir=f0_data_loop_bedpe/data_1st_submission_rep_combined
mkdir -p $outdir
for ii in H3K4me3_Vector H3K4me3_WT H3K4me3_DEL H3K4me3_EIF
do
  for anchor in anchor1 anchor2
  do
    cat ${indir}/${ii}*${anchor}.bed |sort -u  > $outdir/${ii}_${anchor}.bed
  done
done 



# ===========
# batch2 data
# ===========

outdir=f0_data_loop_bedpe/data_202008
mkdir -p $outdir
maps_results_dir=../../f0_looping_data_v2/data_202008_reindex

# anchor1: left anchor, anchor2: right anchor
for bedpe_file in $maps_results_dir/*bedpe
do
    outname=$(basename $bedpe_file .bedpe)
    echo $outname
    tail -n +2 $bedpe_file | awk '{OFS="\t";print$2,$3,$4,$1}' > $outdir/${outname}_anchor1.bed
    tail -n +2 $bedpe_file | awk '{OFS="\t";print$5,$6,$7,$1}' > $outdir/${outname}_anchor2.bed
done




