#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o out_bart3d_dis200k_data_1st_submit.log

# Run program


outdir=bart3d_dis200k_data_1st_submit
mkdir ${outdir}

project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
hicpro_results=${project_dir}/f0_data_process/hichip_hicpro_hichipper/hicpro_data_1st_submit/hicpro/human_hichip/hic_results/matrix/
genomicDistance=200000

# DEL33HICHIP  DEL34HICHIP  EIF3HICHIP  EIF4HICHIP  FL3HICHIP  FL4HICHIP  PCDH3HICHIP  PCDH4HICHIP

treat_arr=(FL  DEL3 EIF)
ctrl_arr=(PCDH FL   DEL3)

for ii in 0 1 2
do
    t1=${hicpro_results}/${treat_arr[ii]}3HICHIP/raw/5000/${treat_arr[ii]}3HICHIP_5000.matrix
    t2=${hicpro_results}/${treat_arr[ii]}4HICHIP/raw/5000/${treat_arr[ii]}4HICHIP_5000.matrix
    c1=${hicpro_results}/${ctrl_arr[ii]}3HICHIP/raw/5000/${ctrl_arr[ii]}3HICHIP_5000.matrix
    c2=${hicpro_results}/${ctrl_arr[ii]}4HICHIP/raw/5000/${ctrl_arr[ii]}4HICHIP_5000.matrix
    bedFileHicpro=${hicpro_results}/${treat_arr[ii]}3HICHIP/raw/5000/${treat_arr[ii]}3HICHIP_5000_abs.bed
    echo $ii
    time bart3d -f hicpro -s hg38 \
           -t $t1,$t2 \
           -c $c1,$c2 \
           --bedFileHicpro $bedFileHicpro \
           --genomicDistance ${genomicDistance} \
           --outdir ${outdir} \
           --outFileName ${treat_arr[ii]}_over_${ctrl_arr[ii]} \
           --coverageNormalization
done


