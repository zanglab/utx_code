#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_out/DEL32.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/dm6"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/dm6/rsem_reference"
mapping_outdir="./processed_results//DEL32"

# ==== STAR MAPPING
STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L001_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L002_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L003_R1_001_val_1.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L004_R1_001_val_1.fq.gz \
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L001_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L002_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L003_R2_001_val_2.fq.gz,\
/nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/trimming_reads/proseq_202101/trim_out/DEL32-PRO_S6_L004_R2_001_val_2.fq.gz \
--outFileNamePrefix ${mapping_outdir}/DEL32 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

# ==== mark duplicates
module load picard

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${mapping_outdir}/DEL32Aligned.sortedByCoord.out.bam \
O=${mapping_outdir}/DEL32Aligned.sortedByCoord.out.markdup.bam \
M=${mapping_outdir}/DEL32Aligned.sortedByCoord.out.markdup.txt \
ASSUME_SORT_ORDER=coordinate 


# ==== RSEM quantification
module load perl

rsem-calculate-expression \
--num-threads 4 \
--fragment-length-max 1000 \
--estimate-rspd \
--no-bam-output \
--paired-end \
--bam ${mapping_outdir}/DEL32Aligned.toTranscriptome.out.bam \
${rsem_reference} ${mapping_outdir}/DEL32_rsem 

