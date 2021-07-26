#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/WT_rep1.out 

salmon quant --gcBias -i /nv/vol190/zanglab/zw5j/data/index/salmon_index/hg38_v92 -l A -r \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FL-RNA_L001.fastq.gz \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FL-RNA_L002.fastq.gz \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FL-RNA_L003.fastq.gz \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FL-RNA_L004.fastq.gz -p 4 -o processed_results//WT_rep1 
