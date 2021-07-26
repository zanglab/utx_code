#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/UTX_FUSIDR_rep1.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/hg38"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/hg38/rsem_reference"
mapping_outdir="processed_results//UTX_FUSIDR_rep1"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FUS-RNA_L001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FUS-RNA_L002.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FUS-RNA_L003.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//FUS-RNA_L004.fastq.gz \
 \
--outFileNamePrefix ${mapping_outdir}/UTX_FUSIDR_rep1 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/UTX_FUSIDR_rep1Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/UTX_FUSIDR_rep1Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/UTX_FUSIDR_rep1Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--strand-specific \
--no-bam-output \
--bam ${mapping_outdir}/UTX_FUSIDR_rep1Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/UTX_FUSIDR_rep1_rsem 

