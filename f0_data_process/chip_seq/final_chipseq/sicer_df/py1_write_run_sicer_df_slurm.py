import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
outdir='sicer_out'

sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202102_H3K27ac_H3K4me1_trim','202011_UTX_trim','202102_UTX_H3K27me3_trim']
# celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
factors= ['UTX','UTXFEB','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3','MLL4']
compr_pairs = [['Vector','WT'],['WT','DEL'],['DEL','EIF']]


for sub_dir in sub_dirs:
    for compr_pair in compr_pairs:
        for factor in factors:
            sicer_outdir='{}/{}_over_{}_{}'.format(outdir,compr_pair[1],compr_pair[0],factor)
            os.makedirs(sicer_outdir,exist_ok=True)
            basename_treatment = '{}_{}'.format(compr_pair[1],factor)
            basename_control = '{}_{}'.format(compr_pair[0],factor)
            bam_control ='{}/f0_data_process/chip_seq/final_chipseq/{}/process_qc_out/{}/{}_treat.bam'.format(project_dir,sub_dir,basename_control,basename_control)
            bam_treatment ='{}/f0_data_process/chip_seq/final_chipseq/{}/process_qc_out/{}/{}_treat.bam'.format(project_dir,sub_dir,basename_treatment,basename_treatment)
            if os.path.isfile(bam_control):
                #print(sub_dir,celltype,factor)
                slurmfile = '{}/{}_over_{}_{}.slurm'.format(slurm_dir,compr_pair[1],compr_pair[0],factor)    
                with open(slurmfile,'w') as slurmout:
                    slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')
            
                    slurmout.write('#SBATCH -o {}/slurm_{}_over_{}_{}.out\n\n'.format(slurm_dir,compr_pair[1],compr_pair[0],factor))
                    slurmout.write('time sicer_df -t \\\n{} \\\n{} \\\n-s hg38 --output_directory {}\n'.format(bam_treatment,bam_control,sicer_outdir))


