#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm2_generate_bw.out

bam_to_bw()
{

filedir=$1 
outdir=$2
prename=$3 
genome_file=$4

bam_file=${filedir}/${prename}/${prename}_treat.bam
prename=${outdir}/${prename}
echo $prename


## == PE process
bedtools bamtobed -bedpe -mate1 -i $bam_file 2> ${prename}.bamtobed.warnings |\
awk '{OFS="\t"} ($9=="+" && $10=="-"){print $1,$2,$6,$7,$8,$9}; ($9=="-"&& $10=="+"){print $1,$5,$3,$7,$8,$9}'|\
sort -k 1,1 -k 2,2g -k 3,3g |grep -v "_" > ${prename}_PE.bed &
wait

## == reverse the strand & get 5 prime 
cat ${prename}_PE.bed |awk '{OFS="\t"} ($6=="+"){print $1,$2,$2+1,"N","0","-"}; ($6=="-"){print $1,$3-1,$3,"N","0","+"}' > ${prename}_5prime_strandRev.bed &
wait

## == Convert to bedGraph
bedtools genomecov -bg -i ${prename}_5prime_strandRev.bed -g $genome_file -strand + > ${prename}_plus.bdg &
bedtools genomecov -bg -i ${prename}_5prime_strandRev.bed -g $genome_file -strand - > ${prename}_minus_abs.bdg &
wait

## == Invert minus strand & normalized by RPM
readCount=`cat ${prename}_PE.bed | grep "" -c` 
cat ${prename}_minus_abs.bdg | awk '{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${prename}_minus.bdg &
wait
cat ${prename}_plus.bdg | awk -v rc=$readCount '{OFS="\t"} {print $1,$2,$3,$4*1000*1000/rc}' > ${prename}_plus.rpm.bdg &
cat ${prename}_minus.bdg | awk -v rc=$readCount '{OFS="\t"} {print $1,$2,$3,$4*1000*1000/rc}' > ${prename}_minus.rpm.bdg &
wait

## == Convert to bigwig
bdg2bw ${prename}_plus.bdg $genome_file &
bdg2bw ${prename}_minus.bdg $genome_file &
bdg2bw ${prename}_plus.rpm.bdg $genome_file &
bdg2bw ${prename}_minus.rpm.bdg $genome_file &
wait

}





## ==== main 

module load samtools

genome_file=/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes
project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/
filedir=${project_dir}/f0_data_process/pro_seq/data_202101_trim/processed_out/

outdir=f2_processed_data_with_dup_sep_rep
mkdir $outdir

## == DEL

for prename in DEL31PRO  DEL32PRO  FL1PRO  FL2PRO  PCDH1PRO  PCDH2PRO
do
bam_to_bw $filedir $outdir $prename $genome_file
done
