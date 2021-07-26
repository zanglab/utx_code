import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
# import matplotlib
# # matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks")
# from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde



    
## ==== main
outdir = 'f3_diff_expr_DCI_RP_WT_over_Vector'
os.makedirs(outdir,exist_ok=True)

project_dir = "/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir = "/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
factors = ['UTX','H3K4me1','H3K4me3','H3K27ac']


df = pd.read_csv('../data/gene_Expression_DCI_RP/fz_data_combined/TPM_DEseq2_DCI_RPKM_RP.csv',index_col=0)
compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']

for batch in ['Batch1','Batch2']:
    for compr_type in compr_types:
        hm_col = 'RP_H3K4me1_{}_log2FC'.format(compr_type)
        dci_col = '{}_H3K4me3_{}_DCI'.format(batch,compr_type)
        expr_col = '{}_log2FoldChange'.format(compr_type)
        hm_thre = np.log2(1.1)
        dci_thre = -np.log10(0.2)
        expr_thre = np.log2(1.1)
        
        df_increased = df[(df[hm_col]>hm_thre)&(df[dci_col]>dci_thre)&(df[expr_col]>expr_thre)]
        df_increased.to_csv('{}/H3K4me1_RP_{}_H3K4me3_DCI_expr_increased_{}.csv'.format(outdir,batch,compr_type))
        df_increased.to_csv('{}/H3K4me1_RP_{}_H3K4me3_DCI_expr_increased_{}.genes.txt'.format(outdir,batch,compr_type),columns=[], header=False)
        
        df_decreased = df[(df[hm_col]<-1*hm_thre)&(df[dci_col]<-1*dci_thre)&(df[expr_col]<-1*expr_thre)]
        df_decreased.to_csv('{}/H3K4me1_RP_{}_H3K4me3_DCI_expr_decreased_{}.csv'.format(outdir,batch,compr_type))
        df_decreased.to_csv('{}/H3K4me1_RP_{}_H3K4me3_DCI_expr_decreased_{}.genes.txt'.format(outdir,batch,compr_type),columns=[], header=False)
    






