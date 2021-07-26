import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect





# ==== main 
    
slurm_dir = 'run_get_RPKM_slurms'
outdir='rpkm_csv'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

bam_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out'
promoter_file='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed'
udhs_file='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_unionDHS_fc4_50merge.bed'


celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']
chip_markers= ['HA','H3K27ac','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']


for chip_marker in chip_markers:
    for celltype in celltypes:
        basename='{}_{}'.format(celltype,chip_marker)
        bamfile='{}/{}/{}_treat.bam'.format(bam_dir,basename,basename)
        
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
            slurmout.write('python get_RPM_RPKM_count_on_regions.py -i {} -t {} -s hg38 -f bam -o {}/{}_on_Promoter.csv\n'.format(promoter_file,bamfile,outdir,basename))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py -i {} -t {} -s hg38 -f bam -m -e 10000 -o {}/{}_on_Promoter_es10kb.csv\n'.format(promoter_file,bamfile,outdir,basename))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py -i {} -t {} -s hg38 -f bam -o {}/{}_on_UDHS.csv\n'.format(udhs_file,bamfile,outdir,basename))
        


    
    
