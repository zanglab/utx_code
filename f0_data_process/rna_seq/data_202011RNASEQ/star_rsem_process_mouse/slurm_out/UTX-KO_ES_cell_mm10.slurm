#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/UTX-KO_ES_cell_mm10.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/mm10"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/mm10/rsem_reference"
mapping_outdir="processed_results//UTX-KO_ES_cell_mm10"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs10_S10_L001_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs10_S10_L002_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs10_S10_L003_R1_001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/202011RNASEQ///bs10_S10_L004_R1_001.fastq.gz \
 \
--outFileNamePrefix ${mapping_outdir}/UTX-KO_ES_cell_mm10 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/UTX-KO_ES_cell_mm10Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/UTX-KO_ES_cell_mm10Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/UTX-KO_ES_cell_mm10Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--no-bam-output \
--bam ${mapping_outdir}/UTX-KO_ES_cell_mm10Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/UTX-KO_ES_cell_mm10_rsem 

