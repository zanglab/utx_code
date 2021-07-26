
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
bart3d_results_dir=${project_dir}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new/bart3d_DCI_rename
promoter_file=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed

outdir=f1_promoter_DCI
mkdir $outdir

for genomic_dis in 200k 500k
do
  for datatype in data_1st_submit data202008
  do 
    subdir=bart3d_dis${genomic_dis}_${datatype}
    outdir_tmp=${outdir}/${subdir}
    mkdir $outdir_tmp
    
    for dci_profile in ${bart3d_results_dir}/${subdir}/*_differential_score_after_coverageNormalization.bed
    do
      prename=$(basename $dci_profile _differential_score_after_coverageNormalization.bed)
      echo $subdir $prename
      time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py \
      -a $promoter_file -b $dci_profile -e1 -1999 -e2 -1 \
      -s hg38 -p $outdir_tmp/${prename}_promoter_DCI.csv
    done    
  
  done
done
