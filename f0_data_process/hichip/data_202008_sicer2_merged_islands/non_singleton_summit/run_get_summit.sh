maps_results_dir="../maps_results/"
out_dir="summit_bedpe"

for hm in K27ACDEL  K27ACEIF  K27ACTPR  K27ACVEC  K27ACWT  K4M3DEL3  K4M3DTPR  K4M3EIF  K4M3PCDH  K4M3WT
# for hm in K27ACDEL
do
    bedpe_file=${maps_results_dir}/${hm}/MAPS_output/${hm}*current/${hm}*sig3Dinteractions.bedpe
    echo $bedpe_file
    cp $bedpe_file ${out_dir}/${hm}.bedpe
    head -n 1 $bedpe_file > ${out_dir}/${hm}_summits.bedpe
    cat $bedpe_file |grep -v Singleton|grep '1$' >> ${out_dir}/${hm}_summits.bedpe
    tail -n +2 ${out_dir}/${hm}_summits.bedpe | awk '{OFS="\t";print$1,$2,$3,NR,$7,$8,$9,$10,$11,$12,$13,$14}' > ${out_dir}/${hm}_summits_anchor1.bed
    tail -n +2 ${out_dir}/${hm}_summits.bedpe | awk '{OFS="\t";print$4,$5,$6,NR,$7,$8,$9,$10,$11,$12,$13,$14}' > ${out_dir}/${hm}_summits_anchor2.bed

done



