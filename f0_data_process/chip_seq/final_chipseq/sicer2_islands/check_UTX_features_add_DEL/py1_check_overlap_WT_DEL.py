import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from matplotlib_venn import venn3,venn2


cellType_colors = {'Vector':'tab:blue',\
                   'WT':'tab:red',\
                   'DEL':'k',\
                   'EIF':'tab:purple',\
                   'TPR':'tab:green',\
                   'MT2':'tab:orange',\
                   'FUS':'tab:gray'}

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

    
wt_specific_df = pd.read_csv('data/WT_NOT_overlap_with_DEL.bed',sep='\t',index_col=3,header=None).loc[:,:3]
del_specific_df = pd.read_csv('data/DEL_NOT_overlap_with_WT.bed',sep='\t',index_col=3,header=None).loc[:,:3]
wt_overlap_del_df = pd.read_csv('data/WT_overlap_with_DEL.bed.uniq',sep='\t',index_col=3,header=None).loc[:,:3]

# ## plot pie chart

a = wt_specific_df.shape[0] 
b = del_specific_df.shape[0] 
c = wt_overlap_del_df.shape[0]
la,lb = 'WT','DEL'


plt.figure(figsize=(2.6,2.6))
out = venn2(subsets=(a,b,c),set_labels=(cellType_labels[la],cellType_labels[lb]),\
            set_colors = (cellType_colors[la],cellType_colors[lb]))    

# for text in out.set_labels:
#     text.set_fontsize(24)
# for text in out.subset_labels:
#     try:
#         text.set_fontsize(22)
#     except:
#         pass
plt.savefig('WT_DEL_overlapped_Venn2.pdf',bbox_inches = 'tight',pad_inches=0.1,transparent=True)
plt.show()
plt.close()

    
    
    
    
    
