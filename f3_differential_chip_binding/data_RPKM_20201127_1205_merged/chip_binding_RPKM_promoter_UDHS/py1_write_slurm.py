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

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
bam_dir='{}/f0_data_process/chip_seq/data_202012/data_20201127CHIP_20201205CHIP_merged/process_qc_out/'.format(project_dir)
promoter_file='{}/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed'.format(project_dir)
udhs_file='{}/f3_differential_chip_binding/data/hg38_unionDHS_fc4_50merge.bed'.format(project_dir)


celltypes = ['PCDH','WT','DEL3','EIF']
chip_markers= ['K4M1','UTX']


for chip_marker in chip_markers:
    for celltype in celltypes:
        basename='{}{}'.format(celltype,chip_marker)
        bamfile='{}/{}/{}_treat.bam'.format(bam_dir,basename,basename)
        print(os.path.isfile(bamfile) ,basename)
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
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_Promoter.csv\n\n'.format(promoter_file,bamfile,outdir,basename))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -m -e 10000 -o {}/{}_on_Promoter_es10kb.csv\n\n'.format(promoter_file,bamfile,outdir,basename))
            slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_UDHS.csv\n'.format(udhs_file,bamfile,outdir,basename))
        


    
    
