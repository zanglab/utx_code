import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import re,bisect



def linear_regression(x,y):
    xmean = np.mean(x)
    ymean = np.mean(y)
    assert len(x)==len(y)
    sum1,sum2=0,0
    x_values,y_values = x.values,y.values
    for i in np.arange(len(x)):
        sum1 += (x_values[i]-xmean)*(y_values[i]-ymean)
        sum2 += np.power((x_values[i]-xmean),2)
    a = sum1/sum2
    b = ymean-a*xmean
    return a,b,xmean,ymean


def return_residual_and_std(df,x,y,a,b):
    # residual of each region/gene
    y_reg = a*x+b
    y_residual = y-y_reg
    miu = np.mean(y_residual)
    std = np.std(y_residual)
    # p-value of each residual
    df['residual']=y_residual
    df['residual_std']=std
    df['residual_pvalue'] = scipy.stats.norm(miu,std).sf(y_residual)
#     df['residual_pvalue'] = pd.concat([df['residual_pvalue'],(1-df['residual_pvalue'])],axis=1).min(axis=1)
    return df
    


# ==== main 
    
outdir = 'f3_highly_variable_chip_binding'
os.makedirs(outdir,exist_ok=True)


chip_markers= ['HA','H3K27ac','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody','Macs2peak','Sicer2peak']
# celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']
count_thre = 0

for chip_marker in chip_markers:
    for rpkm_pattern in rpkm_patterns:
        csv_file='f2_RPKM_collection/{}_{}_TPM.csv'.format(chip_marker,rpkm_pattern)  
        df = pd.read_csv(csv_file,index_col=0)
        df = df[df.sum(axis=1)>count_thre]
        
        x = np.log2(df.mean(axis=1))
        y = np.log2(df.std(axis=1)/df.mean(axis=1))
        # linear fit of log2mean and log2CV
        df['log2mean'] = x
        df['log2CV'] = y
        a,b,xmean,ymean = linear_regression(x,y)
        df = return_residual_and_std(df,x,y,a,b)
        df.to_csv(outdir+os.sep+'{}_{}_TPM_residual.csv'.format(chip_marker,rpkm_pattern))
        # exit()
