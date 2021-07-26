import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn3_plot(hm,cellTypes,indir,outdir,subdir):

    os.makedirs(outdir+os.sep+subdir,exist_ok=True)
    df1 = pd.read_csv('{}/{}/{}_{}.csv'.format(indir,subdir,hm,cellTypes[0]),index_col=0)
    df2 = pd.read_csv('{}/{}/{}_{}.csv'.format(indir,subdir,hm,cellTypes[1]),index_col=0)
    df3 = pd.read_csv('{}/{}/{}_{}.csv'.format(indir,subdir,hm,cellTypes[2]),index_col=0)
    
    a = set(df1.index)
    b = set(df2.index)
    c = set(df3.index)
    set_labels = [cellType_labels[i] for i in cellTypes]
    set_colors = [cellType_colors[i] for i in cellTypes]
    
    plt.figure(figsize=(2.6,2.6))
    out = venn3([a,b,c],set_labels=set_labels,set_colors =set_colors,alpha=0.5)
    plt.title(hm,fontsize=14)
    plt.savefig(outdir+os.sep+subdir+os.sep+'{}.pdf'.format(hm) ,bbox_inches = 'tight',pad_inches=0.1,transparent=True)
    plt.show()
    plt.close()




# ==== main


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


    
indir = '../f2_loop_anchor_annotation/f2_loop_genomic_feature'
outdir = 'f2_loop_venn'
subdirs=['data_1st_submission_rep_combined','data_202008']

# ==== batch1 H3K4me3
subdir=subdirs[0]
hm='H3K4me3'
cellTypes = ['Vector','WT','DEL']
venn3_plot(hm,cellTypes,indir,outdir,subdir)

 
# ==== batch2 H3K4me3
subdir=subdirs[1]
hm='H3K4me3'
venn3_plot(hm,cellTypes,indir,outdir,subdir)

# ==== batch2 H3K27ac
hm='H3K27ac'
venn3_plot(hm,cellTypes,indir,outdir,subdir)


