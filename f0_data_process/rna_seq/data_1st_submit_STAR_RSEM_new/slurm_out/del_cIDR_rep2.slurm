#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/del_cIDR_rep2.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/hg38"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/hg38/rsem_reference"
mapping_outdir="processed_results//del_cIDR_rep2"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//DEL31-RNA_L001.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//DEL31-RNA_L002.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//DEL31-RNA_L003.fastq.gz,\
/nv/vol190/zanglab/bs6yk/NATDATA/RNASEQ//DEL31-RNA_L004.fastq.gz \
 \
--outFileNamePrefix ${mapping_outdir}/del_cIDR_rep2 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/del_cIDR_rep2Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/del_cIDR_rep2Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/del_cIDR_rep2Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--no-bam-output \
--bam ${mapping_outdir}/del_cIDR_rep2Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/del_cIDR_rep2_rsem 

