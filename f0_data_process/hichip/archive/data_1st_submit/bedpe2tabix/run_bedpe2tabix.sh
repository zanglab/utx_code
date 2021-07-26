
maps_results_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_1st_submit/maps_results/


plot_dir=maps_tabix
mkdir $plot_dir

# for mark in K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT
for mark in DEL33HICHIP  DEL34HICHIP  EIF3HICHIP  EIF4HICHIP  FL3HICHIP  FL4HICHIP  PCDH3HICHIP  PCDH4HICHIP
do
    bedpe_file=${maps_results_dir}${mark}/MAPS_output/${mark}*_current/*sig3Dinteractions.bedpe
    echo $mark
    Rscript bedpe2tabix.R -i $bedpe_file -t ${plot_dir}/${mark}
done