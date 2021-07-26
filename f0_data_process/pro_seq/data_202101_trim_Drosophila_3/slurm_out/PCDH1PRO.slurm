#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/PCDH1PRO.out 

module load bwa
module load samtools

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/dm6/dm6.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L001_R1_001_val_1.fq.gz /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L001_R2_001_val_2.fq.gz > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L001.sam
samtools view -bf 0x2 -q 10 -@ 8 ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L001.sam > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/dm6/dm6.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L002_R1_001_val_1.fq.gz /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L002_R2_001_val_2.fq.gz > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L002.sam
samtools view -bf 0x2 -q 10 -@ 8 ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L002.sam > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L002.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/dm6/dm6.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L003_R1_001_val_1.fq.gz /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L003_R2_001_val_2.fq.gz > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L003.sam
samtools view -bf 0x2 -q 10 -@ 8 ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L003.sam > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L003.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/dm6/dm6.fa /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L004_R1_001_val_1.fq.gz /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH1-PRO_S1_L004_R2_001_val_2.fq.gz > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L004.sam
samtools view -bf 0x2 -q 10 -@ 8 ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L004.sam > ./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L004.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g dm -f BAMPE -t \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L001.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L002.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L003.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L004.bam -n ./processed_out/PCDH1PRO/PCDH1PRO

bdg2bw ./processed_out/PCDH1PRO/PCDH1PRO_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/dm6/dm6.chrom.sizes

samtools merge -@ 4 ./processed_out/PCDH1PRO/PCDH1PRO_treat.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L001.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L002.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L003.bam \
./processed_out/PCDH1PRO/PCDH1PRO_PCDH1-PRO_S1_L004.bam

python 2_chipseq_qc_pe.py -i name_list_chipseq_pe.txt -o slurm_out/PCDH1PRO_qc.log -s dm6 -b PCDH1PRO -f1 _R1_001_val_1.fq.gz -f2 _R2_001_val_2.fq.gz

