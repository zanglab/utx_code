#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm3b_generate_bw.out



## ==== main 

module load samtools

genome_file=/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/
filedir=${project_dir}/f0_data_process/pro_seq/data_202101_trim/processed_out/

indir=f3_processed_data_with_dup_sep_rep_MAPQ10
outdir=f3b_merged_plus_minus_tracks
mkdir $outdir


for prename in DEL31PRO  DEL32PRO  FL1PRO  FL2PRO  PCDH1PRO  PCDH2PRO
do
  bedtools unionbedg -i ${indir}/${prename}_plus.rpm.bdg ${indir}/${prename}_minus.rpm.bdg > ${outdir}/${prename}.rpm.union.bdg &
  wait
  cat ${outdir}/${prename}.rpm.union.bdg |awk '{OFS="\t";if ($4>-1*$5) signal=$4;else signal=$5; print$1,$2,$3,signal}' > ${outdir}/${prename}.rpm.bdg &
  wait 
  bdg2bw ${outdir}/${prename}.rpm.bdg $genome_file &
  wait
done
