
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
bart3d_results_dir=${project_dir}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new
utx_dir=${project_dir}//f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding

outdir=f1_UTX_binding_DCI
mkdir $outdir

for genomic_dis in 200k 500k
do
  for datatype in data_1st_submit data202008
  do 
    subdir=bart3d_dis${genomic_dis}_${datatype}
    outdir_tmp=${outdir}/${subdir}
    mkdir $outdir_tmp
    echo $outdir_tmp 
    for dci_profile in ${bart3d_results_dir}/${subdir}/*_differential_score_after_coverageNormalization.bed
    do
      prename=$(basename $dci_profile _differential_score_after_coverageNormalization.bed)
      
      for utx_file_prename in UTX_peaks UTX_islands UTXFEB_peaks UTXFEB_islands
      do
        utx_file=${utx_dir}/${utx_file_prename}.bed
        echo $prename  $utx_file_prename
        
        time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info_revised.py \
        -a $utx_file -b $dci_profile -e2 -1 \
        -s hg38 -p $outdir_tmp/${prename}_${utx_file_prename}_DCI.csv
        
        done
    done    
  done
done
