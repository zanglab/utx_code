import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect


def write_slurm(slurm_dir,fq_dir,basename,fq_file,flag):

    slurmfile = slurm_dir+os.sep+'run_{}_{}.slurm'.format(basename,flag)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}_{}.out\n\n'.format(slurm_dir,basename,flag))
        
        slurmout.write('time split_reads.py -n 10000000 --results_folder check_on_terminal_split_fq/{} {}\n'.format(basename,fq_file))



def main():
    
    slurm_dir = 'run_split_slurms'
    os.makedirs(slurm_dir,exist_ok=True)
    
    fq_dir='/nv/vol190/zanglab/bs6yk/HICHIP/'
    hms = ['DEL33HICHIP',
'DEL34HICHIP',
'EIF3HICHIP',
'EIF4HICHIP',
'FL3HICHIP',
'FL4HICHIP',
'PCDH3HICHIP',
'PCDH4HICHIP']
    for hm in hms:
        files = sorted(glob.glob('{}/{}*.gz'.format(fq_dir,hm)))
        print(hm,files)
        write_slurm(slurm_dir,fq_dir,hm,files[0],'1')
        write_slurm(slurm_dir,fq_dir,hm,files[1],'2')
    


    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
