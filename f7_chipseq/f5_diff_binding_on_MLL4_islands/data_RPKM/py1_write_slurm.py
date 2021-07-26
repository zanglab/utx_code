import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect


def return_chipseq_bam(project_dir,basename):
    # factor in ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
    
    # sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202102_H3K27ac_H3K4me1_trim','202011_UTX_trim']
    sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202102_H3K27ac_H3K4me1_trim','202102_UTX_H3K27me3_trim','202011_UTX_trim']
    for sub_dir in sub_dirs:
        bam_file='{}/f0_data_process/chip_seq/final_chipseq/{}/process_qc_out/{}/{}_treat.bam'.format(project_dir,sub_dir,basename,basename)
        if os.path.isfile(bam_file):
            return bam_file
    return 1
    



# ==== main 
    
slurm_dir = 'run_get_RPKM_slurms'
outdir='rpkm_csv'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
#promoter_file='/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'
#udhs_file='/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
island_factor='MLL4'
islands_file='{}/f0_data_process/chip_seq/final_chipseq/sicer_df/sicer_out/WT_over_Vector_{}/WT_{}_treat-vs-Vector_{}_treat-W200-G600-E1000-union.island'.format(project_dir,island_factor,island_factor,island_factor)

celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
factors= ['UTX','UTXFEB','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3','MLL4']


for factor in factors:
    for celltype in celltypes:
        basename='{}_{}'.format(celltype,factor)
        bamfile = return_chipseq_bam(project_dir,basename)
        #print(basename,bamfile)
        if not os.path.isfile(bamfile):
            continue
        print(basename)
        # ==== write the .slurm file
        slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
        with open(slurmfile,'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
            slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -m -e 500 -o {}/{}_on_{}_islands_es500bp.csv\n\n'.format(islands_file,bamfile,outdir,basename,island_factor))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -m -e 2000 -o {}/{}_on_{}_islands_es2kb.csv\n\n'.format(islands_file,bamfile,outdir,basename,island_factor))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_{}_islands.csv\n\n'.format(islands_file,bamfile,outdir,basename,island_factor))
            # slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -m -e 10000 -o {}/{}_on_Promoter_es10kb.csv\n\n'.format(promoter_file,bamfile,outdir,basename))
            # slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_UDHS.csv\n'.format(udhs_file,bamfile,outdir,basename))
        


    
    
