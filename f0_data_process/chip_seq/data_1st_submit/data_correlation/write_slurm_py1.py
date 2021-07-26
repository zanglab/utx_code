import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

    
slurm_dir = 'slurm_files'
os.makedirs(slurm_dir,exist_ok=True)

outdir = 'RPKM_csv'
os.makedirs(outdir,exist_ok=True)

chipseq_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out'
name_files = pd.read_excel('peak_file_names_correlation.xlsx',index_col=0) 
print(name_files)

for chipseq_name in name_files.index:
    bamfile = '{}/{}/{}_treat.bam'.format(chipseq_dir,chipseq_name,chipseq_name)
#     print(chipseq_name,bamfile)
    slurmfile='{}/{}.slurm'.format(slurm_dir,chipseq_name)
    
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,chipseq_name))
        
        slurmout.write('time python get_RPKM_on_regions_readcount.py -s hg38 -f bam -i hg38_5000_ord.bed \\\n-t {} \\\n-o {}/{}_RPKM_genome_5kb.csv\n\n'.format(bamfile,outdir,chipseq_name))
        slurmout.write('time python get_RPKM_on_regions_readcount.py -s hg38 -f bam -i hg38_4k_promoter_refGene.bed \\\n-t {} \\\n-o {}/{}_RPKM_promoter.csv'.format(bamfile,outdir,chipseq_name))
        

