import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect


def write_slurm(slurm_dir,basename,commandline):

    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
        slurmout.write(commandline)



# ==== main 
    
slurm_dir = 'run_get_RPKM_slurms'
os.makedirs(slurm_dir,exist_ok=True)

bam_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/data_1st_submit/process_qc_out'
peak_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data'
promoter_file='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_4k_promoter_geneID.bed'
gb_file='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f3_differential_chip_binding/data/hg38_minus2kb_and_GB_geneID.bed'

celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']
chip_markers= ['HA','H3K27ac','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']

for chip_marker in chip_markers:
    for celltype in celltypes:
        basename='{}_{}'.format(celltype,chip_marker)
        bamfile='{}/{}/{}_treat.bam'.format(bam_dir,basename,basename)
        peakfile='{}/macs2_peaks_merged/union_{}_peaks.narrowPeak.sorted.merged.id'.format(peak_dir,chip_marker)
        island_file='{}/sicer2_islands_merged/union_{}_treat.scoreisland.sorted.merged.id'.format(peak_dir,chip_marker)
        
        # ==== write the .slurm file
        slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
        with open(slurmfile,'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
            slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
            slurmout.write('python get_RPKM_and_count_on_regions_readcount.py -i {} -t {} -s hg38 -f bam -o rpkm_csv/{}_on_Macs2peak.csv\n'.format(peakfile,bamfile,basename))
            slurmout.write('python get_RPKM_and_count_on_regions_readcount.py -i {} -t {} -s hg38 -f bam -o rpkm_csv/{}_on_Sicer2peak.csv\n'.format(island_file,bamfile,basename))
            slurmout.write('python get_RPKM_and_count_on_regions_readcount.py -i {} -t {} -s hg38 -f bam -o rpkm_csv/{}_on_Promoter.csv\n'.format(promoter_file,bamfile,basename))
            slurmout.write('python get_RPKM_and_count_on_regions_readcount.py -i {} -t {} -s hg38 -f bam -o rpkm_csv/{}_on_GeneBody.csv\n'.format(gb_file,bamfile,basename))
        


    
    
