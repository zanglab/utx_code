
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
# python_file=/nv/vol190/zanglab/zw5j/scripts/Modules/find_overlap_keep_info_NOT_sep_strand.py
bart3d_dir=${project_dir}/f5_hichip/f1_hichip_bart3d/f0_run_bart3d_new
WT_UTX_peak_file=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/WT_UTXFEB_with_Vector_control/WT_UTXFEB_with_Vector_control_q01_peaks.narrowPeak
cat /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/*UTXFEB*/*peaks.narrowPeak > 2021FEB_UTX_all_peaks.bed
all_UTX_peak_file=2021FEB_UTX_all_peaks.bed


for bart_subdir in bart3d_dis200k_data202008 bart3d_dis500k_data202008 bart3d_dis200k_data_1st_submit bart3d_dis500k_data_1st_submit
do
    outdir=overlapped_${bart_subdir}
    mkdir $outdir
    for dci_profile in ${bart3d_dir}/${bart_subdir}/*differential_score*.bed
    do 
      prename=$(basename $dci_profile .bed)
      echo $bart_subdir $prename
      bedtools intersect -a $dci_profile -b $WT_UTX_peak_file -wa -u > ${outdir}/${prename}.WT_UTX_bound.bed
      bedtools intersect -a $dci_profile -b $all_UTX_peak_file -wa -v > ${outdir}/${prename}.all_UTX_unbound.bed
    done
done
     