import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




sub_dirs=['re_1st_submission_H3K4me3_MLL4SC_trim','re_202012_H3K4me2_trim','202104_UTX_H3K4me2','202102_H3K27ac_H3K4me1_trim','202102_UTX_H3K27me3_trim','202011_UTX_trim']
celltypes = ['Vector','WT','DEL','EIF','MT2','TPR']
factors= ['UTX','UTXFEB','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac','MLL4','H3K27me3']

df = pd.DataFrame(columns=['total','total_in_islads','Frip'])
for sub_dir in sub_dirs:
    for celltype in celltypes:
        for factor in factors:
            total_file = 'readCount_csv/{}/{}_{}_on_sicer_islands.csv.total'.format(sub_dir,celltype,factor)
            if os.path.isfile(total_file):
                #print(sub_dir,celltype,factor)
                total_df = pd.read_csv(total_file,sep='\t').rename(index={0:'{}_{}'.format(celltype,factor)})
                df = pd.concat([df,total_df])
df.to_csv('total_reads_in_Islands.csv')

