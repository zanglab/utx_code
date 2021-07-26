#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/EIF_UTXAPR_with_Vector_control.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/EIFUTXAPR_S4_L001_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/EIFUTXAPR_S4_L002_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/EIFUTXAPR_S4_L003_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/EIFUTXAPR_S4_L004_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/PCDHUTXAPR_S1_L001_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/PCDHUTXAPR_S1_L002_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/PCDHUTXAPR_S1_L003_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/2021APRCHIPSEQ/PCDHUTXAPR_S1_L004_R1_001.fastq.gz > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.sam
samtools view -bS -q 1 -@ 8 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.sam > ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.bam -c \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.bam -n ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control

macs2 callpeak --SPMR -B -q 0.1 --keep-dup 1 -g hs -t \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.bam -c \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.bam -n ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_q01

bdg2bw ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

samtools merge -@ 4 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_treat.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_EIFUTXAPR_S4_L004_R1_001.bam

samtools merge -@ 4 ./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_control.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L001_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L002_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L003_R1_001.bam \
./process_qc_out/EIF_UTXAPR_with_Vector_control/EIF_UTXAPR_with_Vector_control_PCDHUTXAPR_S1_L004_R1_001.bam

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/EIF_UTXAPR_with_Vector_control_qc.log -s hg38 -b EIF_UTXAPR_with_Vector_control -f .fastq.gz

