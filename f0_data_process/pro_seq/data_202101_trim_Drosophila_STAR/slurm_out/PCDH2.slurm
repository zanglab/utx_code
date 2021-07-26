#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/PCDH2.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/dm6"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/dm6/rsem_reference"
mapping_outdir="./processed_results//PCDH2"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L001_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L002_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L003_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L004_R1_001_val_1.fq.gz \
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L001_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L002_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L003_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/PCDH2-PRO_S4_L004_R2_001_val_2.fq.gz \
--outFileNamePrefix ${mapping_outdir}/PCDH2 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/PCDH2Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/PCDH2Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/PCDH2Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--no-bam-output \
--paired-end \
--bam ${mapping_outdir}/PCDH2Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/PCDH2_rsem 

