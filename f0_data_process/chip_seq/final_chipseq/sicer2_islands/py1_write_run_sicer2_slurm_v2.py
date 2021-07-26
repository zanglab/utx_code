import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




slurm_dir = 'slurm_files_v2'
os.makedirs(slurm_dir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
outdir='sicer_out'

# sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202102_H3K27ac_H3K4me1_trim','202011_UTX_trim','202102_UTX_H3K27me3_trim']
sub_dirs=['202104_UTX_H3K4me2']
celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
# factors= ['UTX','UTXFEB','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3','MLL4']
factors = ['H3K4me2APR']

for sub_dir in sub_dirs:
    #sicer_outdir='{}/{}'.format(outdir,sub_dir)
    #os.makedirs(sicer_outdir,exist_ok=True)
    for celltype in celltypes:
        for factor in factors:
            basename = '{}_{}'.format(celltype,factor)
            bam_file='{}/f0_data_process/chip_seq/final_chipseq/{}/process_qc_out/{}/{}_treat.bam'.format(project_dir,sub_dir,basename,basename)
            if os.path.isfile(bam_file):
                #print(sub_dir,celltype,factor)
                slurmfile='{}/{}.slurm'.format(slurm_dir,basename)    
                with open(slurmfile,'w') as slurmout:
                    slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')
            
                    slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
                    slurmout.write('time sicer -t {} \\\n-s hg38 --output_directory {}\n'.format(bam_file,outdir))


