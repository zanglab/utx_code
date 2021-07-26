
# data 1st submission
maps_results_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_1st_submit_sicer2_merged_islands_new/maps_results//"
outdir=data_1st_submission
mkdir ${outdir}

hm_array=(DEL33HICHIP  DEL34HICHIP  EIF3HICHIP  EIF4HICHIP  FL3HICHIP  FL4HICHIP  PCDH3HICHIP  PCDH4HICHIP)
rename_array=(H3K4me3_DEL_rep1 H3K4me3_DEL_rep2 H3K4me3_EIF_rep1 H3K4me3_EIF_rep2 H3K4me3_WT_rep1 H3K4me3_WT_rep2 H3K4me3_Vector_rep1 H3K4me3_Vector_rep2)
for ii in {0..7}
do
    hm=${hm_array[ii]}
    rename=${rename_array[ii]}
    bedpe_file=${maps_results_dir}/${hm}/MAPS_output/${hm}*current/${hm}*sig3Dinteractions.bedpe
    echo $bedpe_file
    echo $rename
    cp $bedpe_file ${outdir}/${rename}.bedpe
done



# data 202008_new
maps_results_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_202008_sicer2_merged_islands_with_new_k27/maps_results/"
outdir=data_202008
mkdir ${outdir}

hm_array=(K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT)
rename_array=(H3K27ac_DEL H3K27ac_EIF H3K27ac_TPR H3K27ac_Vector H3K27ac_WT H3K4me3_DEL H3K4me3_TPR H3K4me3_EIF H3K4me3_Vector H3K4me3_WT)
for ii in {0..9}
do
    hm=${hm_array[ii]}
    rename=${rename_array[ii]}
    bedpe_file=${maps_results_dir}/${hm}/MAPS_output/${hm}*current/${hm}*sig3Dinteractions.bedpe
    echo $bedpe_file
    echo $rename
    cp $bedpe_file ${outdir}/${rename}.bedpe
done



