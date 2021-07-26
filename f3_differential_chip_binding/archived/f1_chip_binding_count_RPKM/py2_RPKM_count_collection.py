import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




# ==== main 
    
outdir = 'f2_RPKM_collection'
os.makedirs(outdir,exist_ok=True)


chip_markers= ['HA','H3K27ac','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody','Macs2peak','Sicer2peak']
celltypes = ['Vector','del_cIDR','del_TPR','WT','UTX_eIFIDR']


for chip_marker in chip_markers:
    for rpkm_pattern in rpkm_patterns:
        df_rpkm = pd.DataFrame()
        df_count = pd.DataFrame()
        # collected RPKM for each chip-marker and each rpkm-pattern 
        for celltype in celltypes:
            basename='{}_{}'.format(celltype,chip_marker)
            csv_file='rpkm_csv/{}_on_{}.csv'.format(basename,rpkm_pattern)            
            df = pd.read_csv(csv_file,sep='\t',index_col=0)
            # get the RPKM and read count
            df_rpkm[celltype] = df['RPKM']
            df_count[celltype] = df['ReadCount']
            #print(df);exit(0)
        
        df_rpkm.to_csv(outdir+os.sep+'{}_{}_RPKM.csv'.format(chip_marker,rpkm_pattern))
        df_count.to_csv(outdir+os.sep+'{}_{}_ReadCount.csv'.format(chip_marker,rpkm_pattern))
        # convert RPKM to TPM
        df_tpm = pd.DataFrame.divide(df_rpkm,df_rpkm.sum()/1000000,axis=1)#;print(df_tpm.sum())
        df_tpm.to_csv(outdir+os.sep+'{}_{}_TPM.csv'.format(chip_marker,rpkm_pattern))#;exit()



    
    
