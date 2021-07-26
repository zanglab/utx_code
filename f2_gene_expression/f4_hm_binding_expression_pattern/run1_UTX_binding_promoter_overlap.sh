
tss_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed.tss.sort.merge.bed
# promoter_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
hm_dir=${project_dir}//f7_chipseq/f8_diff_binding_on_k4me1_k27ac_regions/data_hm_binding

outdir=f1_HM_binding_promoter_overlap
mkdir $outdir

      
for hm_file_prename in H3K27ac_peaks H3K4me1_peaks H3K27ac_islands_increased H3K4me1_islands_increased
# for hm_file_prename in UTX_peaks 
do
  hm_file=${hm_dir}/${hm_file_prename}.bed
  for extend_region in 0 2000 10000 50000   
  do  
    time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info_revised.py \
    -a $hm_file -b $tss_file -e1 $extend_region \
    -s hg38 -p $outdir/${hm_file_prename}_tss_es${extend_region}bp.csv
  done       
done
