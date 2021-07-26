#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_files/slurm_EIF_over_DEL_UTXFEB.out

time sicer_df -t \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/EIF_UTXFEB/EIF_UTXFEB_treat.bam \
/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202102_UTX_H3K27me3_trim/process_qc_out/DEL_UTXFEB/DEL_UTXFEB_treat.bam \
-s hg38 --output_directory sicer_out/EIF_over_DEL_UTXFEB
