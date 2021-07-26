#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out_qc/DEL32PRO.out 

module load bwa
module load samtools

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g dm -f BAMPE -t \
./processed_out/DEL32PRO/DEL32PRO_DEL32-PRO_S6_L001.bam \
./processed_out/DEL32PRO/DEL32PRO_DEL32-PRO_S6_L002.bam \
./processed_out/DEL32PRO/DEL32PRO_DEL32-PRO_S6_L003.bam \
./processed_out/DEL32PRO/DEL32PRO_DEL32-PRO_S6_L004.bam -n ./processed_out/DEL32PRO/DEL32PRO

bdg2bw ./processed_out/DEL32PRO/DEL32PRO_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/dm6/dm6.chrom.sizes

python 2_chipseq_qc_pe.py -i name_list_chipseq_pe.txt -o slurm_out_qc/DEL32PRO_qc.log -s dm6 -b DEL32PRO -f1 _R1_001_val_1.fq.gz -f2 _R2_001_val_2.fq.gz

