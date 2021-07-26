import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

    
slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/"
processed_bam_dir='{}/data_1st_submit/process_qc_out'.format(project_dir)

treat_pairs = [['WT_UTX','Vector_UTX'],['del_cIDR_UTX','WT_UTX'],['UTX_eIFIDR_UTX','del_cIDR_UTX']]

for treat_pair in treat_pairs:
    treat1_bam = '{}/{}/{}_treat.bam'.format(processed_bam_dir,treat_pair[0],treat_pair[0])
    treat2_bam = '{}/{}/{}_treat.bam'.format(processed_bam_dir,treat_pair[1],treat_pair[1])
    if os.path.isfile(treat1_bam) and os.path.isfile(treat2_bam):
        slurmfile='{}/{}_vs_{}.slurm'.format(slurm_dir,treat_pair[0],treat_pair[1])
    
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
            
            slurmout.write('#SBATCH -o {}/slurm_{}_{}.out\n\n'.format(slurm_dir,treat_pair[0],treat_pair[1]))
            slurmout.write('processed_bam_dir="{}"\n\n'.format(processed_bam_dir))
        
            slurmout.write('time sicer_df -t {} \\\n{} \\\n-s hg38 --output_directory sicer_results/{}_vs_{}\n'.format(treat1_bam,treat2_bam,treat_pair[0],treat_pair[1]))



