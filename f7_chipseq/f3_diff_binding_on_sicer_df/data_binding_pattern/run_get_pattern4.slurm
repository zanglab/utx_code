#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_get_pattern4.out


vector_k4=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K4me1/Vector_H3K4me1_treat.bam
vector_k27=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/Vector_H3K27ac/Vector_H3K27ac_treat.bam
wt_k4=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K4me1/WT_H3K4me1_treat.bam
wt_k27=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_H3K27ac_H3K4me1_trim/process_qc_out/WT_H3K27ac/WT_H3K27ac_treat.bam
sicer_df_k4_union=data_sicer_union_islands/WT_H3K4me1_treat-vs-Vector_H3K4me1_treat-W200-G600-E1000-union.island.bed
sicer_df_k27_union=data_sicer_union_islands/WT_H3K27ac_treat-vs-Vector_H3K27ac_treat-W200-G600-E1000-union.island.bed


outdir=readCount_csv
mkdir $outdir

# time python get_pattern_near_site_readcount_write_out.py \
# -i $sicer_df_k4_union \
# -t $vector_k4 \
# -s hg38 -f bam -w 2000 -b 200 -m -o $outdir/Vector_H3K4me1_es2kb_bin200.csv &
# 
# time python get_pattern_near_site_readcount_write_out.py \
# -i $sicer_df_k4_union \
# -t $wt_k4 \
# -s hg38 -f bam -w 2000 -b 200 -m -o $outdir/WT_H3K4me1_es2kb_bin200.csv &
# 
# 
# time python get_pattern_near_site_readcount_write_out.py \
# -i $sicer_df_k27_union \
# -t $vector_k27 \
# -s hg38 -f bam -w 2000 -b 200 -m -o $outdir/Vector_H3K27ac_es2kb_bin200.csv &

time python get_pattern_near_site_readcount_write_out.py \
-i $sicer_df_k27_union \
-t $wt_k27 \
-s hg38 -f bam -w 2000 -b 200 -m -o $outdir/WT_H3K27ac_es2kb_bin200.csv 


