
maps_results_dir=summit_bedpe/
plot_dir=maps_tabix
mkdir $plot_dir

for mark in K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT
# for mark in K27ACDEL
do
    bedpe_file=${maps_results_dir}${mark}_summits.bedpe
    echo $mark
    Rscript bedpe2tabix.R -i $bedpe_file -t ${plot_dir}/${mark}
done