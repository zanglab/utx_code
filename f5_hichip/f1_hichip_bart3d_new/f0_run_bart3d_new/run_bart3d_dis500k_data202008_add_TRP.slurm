#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o out_bart3d_dis500k_data202008_add_TPR.log

# Run program


outdir=bart3d_dis500k_data202008
mkdir ${outdir}

project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang
hicpro_results=${project_dir}/f0_data_process/hichip_hicpro_hichipper/hicpro_data_202008/hicpro/human_hichip/hic_results/matrix/
genomicDistance=500000

treat_arr=(K4M3DTPR K27ACTPR)
ctrl_arr=(K4M3WT K27ACWT)

for ii in 0 1
do
    time bart3d -f hicpro -s hg38 \
           -t ${hicpro_results}/${treat_arr[ii]}/raw/5000/${treat_arr[ii]}_5000.matrix \
           -c ${hicpro_results}/${ctrl_arr[ii]}/raw/5000/${ctrl_arr[ii]}_5000.matrix \
           --bedFileHicpro ${hicpro_results}/${treat_arr[ii]}/raw/5000/${treat_arr[ii]}_5000_abs.bed \
           --genomicDistance ${genomicDistance} \
           --outdir ${outdir} \
           --outFileName ${treat_arr[ii]}_over_${ctrl_arr[ii]} \
           --coverageNormalization
done


