#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_proseq2.out



export PATH="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/pro_seq/data_202101_proseq2/seqtk/:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/fastx_bin/bin:$PATH"
module load perl


bwaIndex=/nv/vol190/zanglab/zw5j/data/index/bwa_index/hg38/hg38.fa
chromInfo=/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes 

bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -PE --RNA3=R1_5prime -T myOutput -O myOutput -I test \
--ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC --ADAPT2=TGGAATTCTCGGGTGCCAAGG



