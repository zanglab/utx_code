

# anchor1: left anchor, anchor2: right anchor
maps_results_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_202008_sicer2_merged_islands/maps_results/"

for hm in K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT
# for hm in K27ACDEL
do
    bedpe_file=${maps_results_dir}/${hm}/MAPS_output/${hm}*current/${hm}*sig3Dinteractions.bedpe
    echo $bedpe_file
    cp $bedpe_file ${hm}.bedpe
    tail -n +2 ${hm}.bedpe | awk '{OFS="\t";print$1,$2,$3,NR,$7,$8,$9,$10,$11,$12,$13,$14}' > ${hm}_anchor1.bed
    tail -n +2 ${hm}.bedpe | awk '{OFS="\t";print$4,$5,$6,NR,$7,$8,$9,$10,$11,$12,$13,$14}' > ${hm}_anchor2.bed

done

