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
bam_dir='{}/f0_data_process/chip_seq/final_chipseq/re_202012_H3K4me2_trim/process_qc_out/'.format(project_dir)
promoter_file='/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'
udhs_file='/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'


celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
chip_markers= ['H3K4me2']


for chip_marker in chip_markers:
    for celltype in celltypes:
        basename='{}_{}'.format(celltype,chip_marker)
        bamfile='{}/{}/{}_treat.bam'.format(bam_dir,basename,basename)
        if not os.path.isfile(bamfile):
            continue
        print(basename)
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
#             slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_Promoter.csv\n\n'.format(promoter_file,bamfile,outdir,basename))
#             slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_UDHS.csv\n'.format(udhs_file,bamfile,outdir,basename))
            slurmout.write('python get_pattern_near_site_readcount_write_out.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -w 2000 -b 200 -m -o {}/{}_es2kb_bin200_on_Promoter.csv\n\n'.format(promoter_file,bamfile,outdir,basename))
            slurmout.write('python get_pattern_near_site_readcount_write_out.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -w 2000 -b 200 -m -o {}/{}_es2kb_bin200_on_UDHS.csv\n\n'.format(udhs_file,bamfile,outdir,basename))


    
    
