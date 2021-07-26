import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde



def return_differential_binding(factor,indir):
    if factor=='UTX':
        genomic_region='es500bp'
        norm_col='total'
    else:
        genomic_region='es2kb'
        norm_col='total_in_islads'
        
    csv_file='{}/{}_{}_differential_NormReadCount_by_{}.csv'.format(indir,factor,genomic_region,norm_col)
    df_factor = pd.read_csv(csv_file,index_col=0)
    return df_factor



project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
genomic_regions = ['es500bp','es2kb',]
norm_cols = ['total','total_in_islads']


indir='f1_differential_NormReadCount'
outdir='f2_combined_data'
os.makedirs(outdir,exist_ok=True)

df = pd.DataFrame()
for factor in factors:
    df_factor = return_differential_binding(factor,indir)
    df = pd.concat([df,df_factor],axis=1)
# df.to_csv(outdir+os.sep+'combined_differential_NormReadCount.csv')
df.round(8).to_csv(outdir+os.sep+'combined_differential_NormReadCount.csv')


