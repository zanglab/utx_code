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

    
# wt_df = pd.read_csv('data/WT_peaks.bed',sep='\t',index_col=3,header=None).loc[:,:3]
# del_df = pd.read_csv('data/DEL_peaks.bed',sep='\t',index_col=3,header=None).loc[:,:3]
# overlapped_df = pd.read_csv('data/DEL_overlap_with_WT.bed',sep='\t',index_col=3,header=None).loc[:,:3]

# ## plot pie chart of DEL

a,b = 5414, 8359-5414
labels = 'overlap with\nactivate marks', 'Not overlap with \nactivate marks', 
title = 'DEL specific'
fig = plt.figure(figsize=(2.,2.))
sizes = [a,b]
explode = (0,0)
plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90,counterclock=False,colors=('salmon','lightgrey'))
plt.axes().axis('equal') 
plt.title('{}'.format(title))
plt.savefig('{}_overlap_with_HM.pdf'.format(title),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
plt.show()
plt.close()
    

# WT
a,b = 10808, 15647-10808
labels = 'overlap with\nactivate marks', 'Not overlap with \nactivate marks', 
title = 'WT'
fig = plt.figure(figsize=(2.,2.))
sizes = [a,b]
explode = (0,0)
plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90,counterclock=False,colors=('salmon','lightgrey'))
plt.axes().axis('equal') 
plt.title('{}'.format(title))
plt.savefig('{}_overlap_with_HM.pdf'.format(title),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
plt.show()
plt.close()

    
    
    
    
    
