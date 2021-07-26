import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect





# ==== main 
    
slurm_dir = 'run_get_RPKM_at_PROseq_union_peak'
outdir='f5_proseq_union_peak_RPKM'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
bed_dir='{}/f0_data_process/pro_seq/data_202101_trim_bigwig/f3_processed_data_with_dup_sep_rep_MAPQ10'.format(project_dir)

# data_dir='/nv/vol190/zanglab/zw5j/data/'
# promoter_file='{}/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'.format(data_dir)
# gb_file='{}/geneID_annotation/hg38/hg38_plus500bp_and_GB_geneID.bed'.format(data_dir)
# udhs_file='{}/unionDHS/hg38_unionDHS_fc4_50merge.bed'.format(data_dir)
peak_file='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f6_proseq/data_0x2_MAPQ10/data/union_PROseq_peak/proseq_all_peaks.sorted.merged.id.bed'


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
        slurmout.write('python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -o {}/{}_on_union_Proseq_peak.csv\n\n'.format(peak_file,bedfile,outdir,basename))
        slurmout.write('python ../data_modules_revised/get_RPM_RPKM_count_on_regions.py \\\n-i {} \\\n-t {} \\\n-s hg38 -f bed -g 0 -m -e 2000 -o {}/{}_on_union_Proseq_peak_es2kb.csv\n\n'.format(peak_file,bedfile,outdir,basename))


    
    
