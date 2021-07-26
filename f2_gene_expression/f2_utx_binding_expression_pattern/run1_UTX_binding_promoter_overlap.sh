
tss_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed.tss.sort.merge.bed
# promoter_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
utx_dir=${project_dir}//f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding

outdir=f1_UTX_binding_promoter_overlap
mkdir $outdir

      
for utx_file_prename in UTX_peaks UTX_islands UTXFEB_peaks UTXFEB_islands
# for utx_file_prename in UTX_peaks 
do
  utx_file=${utx_dir}/${utx_file_prename}.bed
  for extend_region in 0 2000 10000 50000   
  do  
    time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info_revised.py \
    -a $utx_file -b $tss_file -e1 $extend_region \
    -s hg38 -p $outdir/${utx_file_prename}_tss_es${extend_region}bp.csv
  done       
done
