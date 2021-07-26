#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run1_overlap.log



python_script="find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py"
promoter_file="/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed"
bart3d_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f5_hichip/f1_hichip_bart3d/f0_run_bart3d/

treat_arr=(K4M3WT  K4M3DEL3 K4M3EIF  K27ACWT  K27ACDEL K27ACEIF)
ctrl_arr=(K4M3PCDH K4M3WT   K4M3DEL3 K27ACVEC K27ACWT  K27ACDEL)

for dis_para in dis200k dis500k
do
    outdir=promoter_DCI_${dis_para}
    mkdir ${outdir}
    for ii in 0 1 2 3 4 5
    do
        prename=${treat_arr[ii]}_over_${ctrl_arr[ii]}
        bart3d_result=${bart3d_dir}/bart3d_${dis_para}/${prename}_differential_score.bed
        echo $dis_para $prename
        time python $python_script -a $promoter_file -b $bart3d_result -s hg38 -e2 -1 -p ${outdir}/${prename}_DCI_on_promoter.bed
    done
done




