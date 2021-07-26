import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect



# ==== main 
    
slurm_dir = 'run_get_readCount_slurms'
outdir='readCount_csv'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
proseq_dir='{}/f0_data_process/pro_seq/data_202101_trim_bigwig/f4_strand_specific_bedfile'.format(project_dir)
utx_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding'.format(project_dir)

# celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
# factors= ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
utx_prenames = ['UTX_peaks','UTX_islands','UTXFEB_peaks','UTXFEB_islands']

proseq_prenames= {'DEL31PRO':'DEL_rep1',
                  'DEL32PRO':'DEL_rep2',
                  'FL1PRO':'WT_rep1',
                  'FL2PRO':'WT_rep2',
                  'PCDH1PRO':'Vector_rep1',
                  'PCDH2PRO':'Vector_rep2'}


for utx_prename in utx_prenames:
    peak_file='{}/{}.bed'.format(utx_dir,utx_prename)
    for proseq_key in proseq_prenames.keys():
        basename = proseq_prenames[proseq_key]
        proseq_file_plus = '{}/{}_5prime_plus.bed'.format(proseq_dir,proseq_key)
        proseq_file_minus = '{}/{}_5prime_minus.bed'.format(proseq_dir,proseq_key)
        slurmfile = slurm_dir+os.sep+'run_PROseq_{}_on_{}.slurm'.format(basename,utx_prename)
        with open(slurmfile,'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
            slurmout.write('#SBATCH -o {}/slurm_PROseq_{}_on_{}.out\n\n'.format(slurm_dir,basename,utx_prename))
            slurmout.write('python get_pattern_near_site_readcount_write_out_revised.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o {}/PROseq_{}_plus_on_{}_es2kb_bin200.csv\n\n'.format(peak_file,proseq_file_plus,outdir,basename,utx_prename))
            slurmout.write('python get_pattern_near_site_readcount_write_out_revised.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o {}/PROseq_{}_minus_on_{}_es2kb_bin200.csv\n\n'.format(peak_file,proseq_file_minus,outdir,basename,utx_prename))

    
    
