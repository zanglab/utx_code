#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/UTX_eIFIDR_H3K4me1_rep.out 

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/NATDATA/CHIPSEQ//EIFK4M12_L001.fastq.gz > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L001.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L001.sam > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L001.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/NATDATA/CHIPSEQ//EIFK4M12_L002.fastq.gz > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L002.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L002.sam > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L002.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/NATDATA/CHIPSEQ//EIFK4M12_L003.fastq.gz > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L003.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L003.sam > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L003.bam

bwa mem -t 8 /nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/bs6yk/NATDATA/CHIPSEQ//EIFK4M12_L004.fastq.gz > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L004.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L004.sam > /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L004.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L001.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L002.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L003.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_EIFK4M12_L004.bam -n /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep

bdg2bw /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out/UTX_eIFIDR_H3K4me1_rep/UTX_eIFIDR_H3K4me1_rep_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/UTX_eIFIDR_H3K4me1_rep_qc.log -s hg38 -b UTX_eIFIDR_H3K4me1_rep -f .fastq.gz

