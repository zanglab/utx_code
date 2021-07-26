
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
bart3d_results_dir=${project_dir}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new
hm_dir=${project_dir}//f7_chipseq/f8_diff_binding_on_k4me1_k27ac_regions/data_hm_binding

outdir=f1_HM_binding_DCI
mkdir $outdir

for genomic_dis in 200k 500k
do
  for datatype in data_1st_submit data202008
  do 
    subdir=bart3d_dis${genomic_dis}_${datatype}
    outdir_tmp=${outdir}/${subdir}
    mkdir $outdir_tmp
    echo $outdir_tmp 
    for dci_profile in ${bart3d_results_dir}/${subdir}/*.bed
    do
      prename=$(basename $dci_profile .bed)
      for hm_file_prename in H3K27ac_peaks H3K4me1_peaks H3K27ac_islands_increased H3K4me1_islands_increased
      do
        hm_file=${hm_dir}/${hm_file_prename}.bed
        echo $prename  $hm_file_prename
        time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info_revised.py \
        -a $hm_file -b $dci_profile -e2 -1 \
        -s hg38 -p $outdir_tmp/${prename}_${hm_file_prename}_DCI.csv
        
      done
    done    
  done
done
