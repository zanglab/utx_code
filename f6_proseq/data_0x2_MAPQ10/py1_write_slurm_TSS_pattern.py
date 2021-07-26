import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect





# ==== main 
    
slurm_dir = 'run_get_pattern_slurms'
outdir='f1_TSS_pattern'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
bed_dir='{}/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10'.format(project_dir)

data_dir='/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38'
promoter_file='{}/hg38_4k_promoter_geneID.bed'.format(data_dir)
# udhs_file='{}/f3_differential_chip_binding/data/hg38_unionDHS_fc4_50merge.bed'.format(project_dir)


basenames= ['DEL31PRO','DEL32PRO','FL1PRO','FL2PRO','PCDH1PRO','PCDH2PRO']


for basename in basenames:
    bedfile='{}/{}_5prime_strandRev.bed'.format(bed_dir,basename)
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
        slurmout.write('python ../data_modules_revised/get_pattern_near_site_RPKM.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o {}/{}_TSS_es2kb_bin200_RPKM.csv\n\n'.format(promoter_file,bedfile,outdir,basename))
        slurmout.write('python ../data_modules_revised/get_pattern_near_site_raw_count.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -w 2000 -b 200 -m -o {}/{}_TSS_es2kb_bin200_rawCount.csv\n\n'.format(promoter_file,bedfile,outdir,basename))
        # slurmout.write('python get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bam -o {}/{}_on_UDHS.csv\n'.format(udhs_file,bamfile,outdir,basename))
        


    
    
