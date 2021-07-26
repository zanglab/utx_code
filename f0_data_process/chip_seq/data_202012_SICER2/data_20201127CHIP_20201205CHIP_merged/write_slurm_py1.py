import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

    
slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/"
processed_bam_dir='{}/data_202012/data_20201127CHIP_20201205CHIP_merged/process_qc_out'.format(project_dir)


for sub_folder in glob.glob('{}/*'.format(processed_bam_dir)):
    if os.path.isdir(sub_folder):
        pre_name = os.path.basename(sub_folder)
        bamfile = glob.glob('{}/{}*_treat.bam'.format(sub_folder,pre_name))[0]
        bamfile = os.path.basename(bamfile)
        print(pre_name)
        
        slurmfile='{}/{}.slurm'.format(slurm_dir,pre_name)
    
        with open(slurmfile,'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
            slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,pre_name))
            slurmout.write('processed_bam_dir="{}"\n\n'.format(sub_folder))
        
            if re.search('Vector_control',pre_name):
                ctrl_bamfile = '{}_control.bam'.format(pre_name)
                slurmout.write('time sicer -t ${{processed_bam_dir}}/{} \\\n-c ${{processed_bam_dir}}/{} \\\n-s hg38 --output_directory sicer_results/{}\n'.format(bamfile,ctrl_bamfile,pre_name))
            else:
                slurmout.write('time sicer -t ${{processed_bam_dir}}/{} \\\n-s hg38 --output_directory sicer_results/{}\n'.format(bamfile,pre_name))
