#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_generate_bw.out2

bam_to_bw()
{

filedir=$1 
rep1=$2
rep2=$3 
outdir=$4 
prename=$5 
genome_file=$6

rep1_bam=${filedir}/${rep1}/${rep1}_treat.bam
rep2_bam=${filedir}/${rep2}/${rep2}_treat.bam
prename=${outdir}/${prename}


## == Merge and filter 
samtools merge -@ 4 ${prename}_raw.bam $rep1_bam $rep2_bam &
wait
samtools view -H ${prename}_raw.bam > ${prename}_filtered.sam &
wait
samtools view -f 0x2 ${prename}_raw.bam |awk 'NR%2==1{mapq=$5;forward=$0} NR%2==0{if($5>=10 && mapq>=10 && substr($3,1,3)=="chr" && $7=="=" ) print forward"\n"$0}' >>${prename}_filtered.sam & 
wait
samtools view -bS ${prename}_filtered.sam > ${prename}.bam &
wait
# rm ${prename}_filtered.sam


## == PE process
bedtools bamtobed -bedpe -mate1 -i ${prename}.bam 2> ${prename}.bamtobed.warnings |\
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
filedir=${project_dir}/f0_data_process/pro_seq/data_202101/processed_out/

outdir=f2_processed_data_with_dup
mkdir $outdir

## == DEL
rep1=DEL31PRO  
rep2=DEL32PRO
prename=del_PRO
bam_to_bw $filedir $rep1 $rep2 $outdir $prename $genome_file


## == WT
rep1=FL1PRO  
rep2=FL2PRO  
prename=wt_PRO
bam_to_bw $filedir $rep1 $rep2 $outdir $prename $genome_file


## == Vector
rep1=PCDH1PRO  
rep2=PCDH2PRO
prename=vector_PRO
bam_to_bw $filedir $rep1 $rep2 $outdir $prename $genome_file

