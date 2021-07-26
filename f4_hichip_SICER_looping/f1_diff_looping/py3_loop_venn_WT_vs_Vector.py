import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn2_plot(a,b,la,lb,figname,factor,color1='tab:red',color2='tab:blue'):
    # fisher exact test
    # total=19826
    # s,p = stats.fisher_exact([[len(a.intersection(b)),len(a)],[len(b),total-len(b)]])
    # print(compr_type,deg_type,len(a),len(b))
    # print('odds ratio={:.2f}, pvalue={:.2e}'.format(s,p))
    
    plt.figure(figsize=(3,3))
    out = venn2([a,b],set_labels=(la,lb),set_colors=(color1,color2), alpha=0.5)
    
    for text in out.set_labels:
        text.set_fontsize(16)
    for text in out.subset_labels:
        try:
            text.set_fontsize(16)
        except:
            pass
    plt.title(factor)
    # if p<0.05:
        # plt.text(x=.2,y=-.15,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches=0.1,transparent=True)
    plt.show()
    plt.close()



# ==== main


outdir = 'f3_loop_venn'
os.makedirs(outdir,exist_ok=True)
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

subdirs=['f1_diff_looping','f2_diff_looping_high_thre']
factors = ['H3K4me3','H3K27ac']

count_df = pd.DataFrame()
for subdir in subdirs[:]:
    for factor in factors[:]:
        prename='{}_{}'.format(subdir,factor)
        csv_file='{}/{}.csv'.format(subdir,factor)
        df = pd.read_csv(csv_file)
        wt_index = df[df['cell_type'].str.contains('WT')].index
        vector_index = df[df['cell_type'].str.contains('VEC')].index
        shared = wt_index.intersection(vector_index)
        count_df.loc[prename,'WT specific loops'] = len(wt_index.difference(shared))
        count_df.loc[prename,'Vector specific loops'] = len(vector_index.difference(shared))
        count_df.loc[prename,'shared loops'] = len(shared)
        
        a = set(vector_index)
        b = set(wt_index)
        la = 'Vector'
        lb = 'WT'
        figname = outdir+os.sep+'{}_venn.png'.format(prename)
        venn2_plot(a,b,la,lb,figname,factor)

count_df.to_csv(outdir+os.sep+'counting_summary.csv')       

