
deg_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"

for hm in H3K27ac H3K4me3
do
    loop_file=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f4_hichip_SICER_looping/f1_diff_looping/f1_diff_looping/${hm}.csv
    for deg in up dn
    do
        deg_file=${deg_dir}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f2_deg/treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_${deg}genes.txt
        # output file
        outfile=f1_loop_anchor_with_deg/${hm}_${deg}.csv
        head -n 1 $loop_file > ${outfile}
        grep -f $deg_file -w $loop_file >> ${outfile}
    done
done
