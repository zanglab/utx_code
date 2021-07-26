#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o bart2_peak_chipseq_1st_submit.out


file_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out"

outdir=bart2_peak_chipseq_1st_submit
mkdir $outdir

for infile in ${file_dir}/*/*peaks.narrow*
do 
  echo $infile
  bart2 region -i $infile -s hg38 -c 9 --outdir $outdir
done




