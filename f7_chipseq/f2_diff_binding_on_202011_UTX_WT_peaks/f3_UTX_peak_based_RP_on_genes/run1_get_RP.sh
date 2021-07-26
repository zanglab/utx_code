peak_file=/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202011_UTX_trim/process_qc_out/WT_UTX_with_Vector_control/WT_UTX_with_Vector_control_peaks.narrowPeak 
ucsc_file=/nv/vol190/zanglab/zw5j/data/ucsc/hg38_unique_geneSymbol.ucsc

python2 get-regulatory-potential-on-genes_peak_level.py \
-s hg38 -b $peak_file -g $ucsc_file -o gene_rp_from_UTX_WT_peak.txt