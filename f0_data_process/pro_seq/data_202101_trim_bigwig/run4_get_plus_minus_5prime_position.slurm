#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm4.out



## ==== main 


indir=f3_processed_data_with_dup_sep_rep_MAPQ10
outdir=f4_strand_specific_bedfile
mkdir $outdir


for prename in DEL31PRO  DEL32PRO  FL1PRO  FL2PRO  PCDH1PRO  PCDH2PRO
do
  cat ${indir}/${prename}_5prime_strandRev.bed |grep -w + > ${outdir}/${prename}_5prime_plus.bed
  cat ${indir}/${prename}_5prime_strandRev.bed |grep -w - > ${outdir}/${prename}_5prime_minus.bed
  
  echo $prename
  wc -l ${indir}/${prename}_5prime_strandRev.bed
  wc -l ${outdir}/${prename}_5prime_plus.bed
  wc -l ${outdir}/${prename}_5prime_minus.bed
done
