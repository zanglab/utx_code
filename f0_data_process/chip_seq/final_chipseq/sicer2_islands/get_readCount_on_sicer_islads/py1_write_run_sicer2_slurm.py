import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
outdir='readCount_csv'

sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202102_H3K27ac_H3K4me1_trim','202011_UTX_trim','202102_UTX_H3K27me3_trim']
celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
factors= ['UTX','UTXFEB','H3K4me1','H3K4me2','H3K4me3','H3K27ac','MLL4','H3K27me3']

for sub_dir in sub_dirs:
    sub_outdir='{}/{}'.format(outdir,sub_dir)
    os.makedirs(sub_outdir,exist_ok=True)
    for celltype in celltypes:
        for factor in factors:
            basename = '{}_{}'.format(celltype,factor)
            bam_file='{}/f0_data_process/chip_seq/final_chipseq/{}/process_qc_out/{}/{}_treat.bam'.format(project_dir,sub_dir,basename,basename)
            island_file='{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/merged_islands/{}.merged.bed'.format(project_dir,factor)
            if os.path.isfile(bam_file):
                #print(sub_dir,celltype,factor)
                slurmfile='{}/{}_{}.slurm'.format(slurm_dir,sub_dir,basename)    
                with open(slurmfile,'w') as slurmout:
                    slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')
            
                    slurmout.write('#SBATCH -o {}/slurm_{}_{}.out\n\n'.format(slurm_dir,sub_dir,basename))
                    #slurmout.write('time sicer -t {} \\\n-s hg38 --output_directory {}\n'.format(bam_file,sicer_outdir))
                    slurmout.write('python get_readCount_on_region_from_BamBed.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_sicer_islands.csv\n\n'.format(island_file,bam_file,sub_outdir,basename))
            


