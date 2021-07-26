import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

    
slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

processed_bam_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out'

name_files = pd.read_excel('peak_bam_file_names.xlsx',index_col=0) 
print(name_files)

for pre_name in name_files.index:
    bamfile = '{}/{}_treat.bam'.format(pre_name,pre_name)
    if re.search('Vector_ctrl',pre_name):
        ctrl_bamfile = '{}/{}_control.bam'.format(pre_name,pre_name)
    
    slurmfile='{}/{}.slurm'.format(slurm_dir,pre_name)
    
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,pre_name))
        slurmout.write('processed_bam_dir="{}"\n\n'.format(processed_bam_dir))
        
        if re.search('Vector_ctrl',pre_name):
            slurmout.write('time sicer -t ${{processed_bam_dir}}/{} \\\n-c ${{processed_bam_dir}}/{} \\\n-s hg38 --output_directory sicer_results/{}\n'.format(bamfile,ctrl_bamfile,pre_name))
        else:
            slurmout.write('time sicer -t ${{processed_bam_dir}}/{} \\\n-s hg38 --output_directory sicer_results/{}\n'.format(bamfile,pre_name))
