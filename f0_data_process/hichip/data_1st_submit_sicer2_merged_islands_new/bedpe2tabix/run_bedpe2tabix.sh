
maps_results_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_1st_submit_sicer2_merged_islands_new/maps_results/
plot_dir=maps_tabix
mkdir $plot_dir

for mark in DEL33HICHIP  DEL34HICHIP  EIF3HICHIP  EIF4HICHIP  FL3HICHIP  FL4HICHIP  PCDH3HICHIP  PCDH4HICHIP
# for mark in K27ACDEL
do
    bedpe_file=${maps_results_dir}${mark}/MAPS_output/${mark}*_current/*sig3Dinteractions.bedpe
    echo $mark
    Rscript bedpe2tabix.R -i $bedpe_file -t ${plot_dir}/${mark}
done