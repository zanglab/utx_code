# README

- This repository contains all source code for analyzing the data and generating the figures in the paper:
  **UTX condensation underlies its tumor suppressive activity**
- A README file is attached to each folder/subfolder.



# Notes of data analysis


### f0_data_process
- QC and process the raw ChIP-seq, RNA-seq and HiChIP datasets, etc. 


### f2_gene_expression
- identify DEG 
- PC analysis


### f3_differential_chip_binding
- differential ChIP-seq binding union peaks/islands, gene promoter (TSS+/-2kb, TSS+/-10kb) and gene body
- then check the binding patterns of DEG between Vector and WT


### f4_hichip_SICER_looping
- HiChIP data processed using MAPS
- check the length of loops and genomic positions of loop anchors
- check the association between genes and differential loops


### f4_hichip_MAPS_new
- HiChIP data processed using MAPS
- using the final selected chip-seq data to call SICER islands
- and use the merged islands as input for MAPS


### f5_hichip
- additional hichip analysis
- HiC-Pro + BART3D


### f6_proseq
- analysis of proseq data

 
### f7_chipseq
- analysis of final selected QC passed chipseq data
- raw & trimmed data /nv/vol190/zanglab/zw5j/projects_data/UTX_HaoJiang/final_chipseq_merge_trim
- processed data /nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq


### f8_integrative_analysis
- integrate expression data, chip-seq binding



<!-- # Notes of data visualization -->




