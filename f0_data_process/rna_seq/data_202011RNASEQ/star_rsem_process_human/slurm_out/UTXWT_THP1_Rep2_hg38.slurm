#!/bin/bash
#SBATCH -n 1 
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/UTXWT_THP1_Rep2_hg38.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/hg38"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/hg38/rsem_reference"
mapping_outdir="processed_results//UTXWT_THP1_Rep2_hg38"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs17_S17_L001_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs17_S17_L002_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs17_S17_L003_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs17_S17_L004_R1_001.fastq.gz \
 \
--outFileNamePrefix ${mapping_outdir}/UTXWT_THP1_Rep2_hg38 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/UTXWT_THP1_Rep2_hg38Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/UTXWT_THP1_Rep2_hg38Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/UTXWT_THP1_Rep2_hg38Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--no-bam-output \
--bam ${mapping_outdir}/UTXWT_THP1_Rep2_hg38Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/UTXWT_THP1_Rep2_hg38_rsem 

