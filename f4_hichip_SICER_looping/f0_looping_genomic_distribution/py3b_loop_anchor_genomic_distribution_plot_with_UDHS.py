import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

    


# plot the figs

def bar_plot(hms,ticklabels,hm_mark,outdir,color):
    plt.figure(figsize=(3,2.6))
    labels = ["P-P","P-E","E-E","Desert"] 
    position = 1
    for hm in hms:
        a = sum_df.loc[labels[0],'{}%'.format(hm)]  
        b = sum_df.loc[labels[1],'{}%'.format(hm)]
        c = sum_df.loc[labels[2],'{}%'.format(hm)]
        d = sum_df.loc['Total','{}%'.format(hm)] -a-b-c
    
        g3 = plt.bar(position,100*d,bottom=0,width = .68, lw=0,color = color,alpha=1,label = labels[3])
        g2 = plt.bar(position,100*c,bottom=100*d,width = .68, lw=0,color = color,alpha=.6,label = labels[2])
        g1 = plt.bar(position,100*b,bottom=100*(d+c),width = .68, lw=0,color = color,alpha=.4,label = labels[1])
        g0 = plt.bar(position,100*a,bottom=100*(d+c+b),width = .68, lw=0,color = color,alpha=.2,label = labels[0])
        position+=1
        
    plt.ylabel('Genomic interactions (%)',fontsize=16)
    plt.ylim([0,100])
    plt.xlim([0.3,6.7])
    plt.title('{}\n'.format(hm_mark),fontsize=16,va='center')
    plt.legend([g0,g1,g2,g3],labels,loc=1,bbox_to_anchor=(1.32,1.05),fontsize=14,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,frameon=False)
    plt.axes().set_xticks([1.2,2.2,3.2,4.2,5.2]) 
    sns.despine(offset=0, trim=True)
    plt.axes().spines['bottom'].set_visible(False)
    plt.axes().set_xticklabels(ticklabels,rotation=45,ha='right',fontsize=15) 
    plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
#     plt.show()
#     plt.savefig('{}/{}_genomic_distribution.pdf'.format(outdir,hm_mark),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.savefig('{}/{}_genomic_distribution.tiff'.format(outdir,hm_mark),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)



# ==== main()

indir = 'f2_loop_genomic_feature_new_with_UDHS'
outdir = 'f3_loop_genomic_distribution_new_with_UDHS'
os.makedirs(outdir,exist_ok=True)

sum_df = pd.read_csv(outdir+os.sep+'loop_distribution_summary.csv',index_col=0)

# plot the figs
ticklabels = ['Vector','WT','$\Delta$cIDR','UTX_eIF$_{IDR}$','$\Delta$TPR']
hms = ['K27ACVEC','K27ACWT','K27ACDEL','K27ACEIF','K27ACTPR',]    
color = 'navy'
bar_plot(hms,ticklabels,'H3K27ac',outdir,color)   

hms = ['K4M3PCDH','K4M3WT','K4M3DEL3','K4M3EIF','K4M3DTPR',]
bar_plot(hms,ticklabels,'H3K4me3',outdir,color)   
    
