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
    labels = ["P-P","P-GB","P-IG","GB-GB","GB-IG","IG-IG"][::-1]
    position = 1
    for hm in hms:
        a = sum_df.loc[labels[0],'{}%'.format(hm)] # a at the bottom
        b = sum_df.loc[labels[1],'{}%'.format(hm)]
        c = sum_df.loc[labels[2],'{}%'.format(hm)]
        d = sum_df.loc[labels[3],'{}%'.format(hm)]
        e = sum_df.loc[labels[4],'{}%'.format(hm)]
        f = sum_df.loc[labels[5],'{}%'.format(hm)]
    
        g0 = plt.bar(position,100*a,bottom=0,width = .68, lw=0,color = color,alpha=1,label = labels[0])
        g1 = plt.bar(position,100*b,bottom=100*a,width = .68, lw=0,color = color,alpha=.85,label = labels[1])
        g2 = plt.bar(position,100*c,bottom=100*(a+b),width = .68, lw=0,color = color,alpha=.7,label = labels[2])
        g3 = plt.bar(position,100*d,bottom=100*(a+b+c),width = .68, lw=0,color = color,alpha=.5,label = labels[3])
        g4 = plt.bar(position,100*e,bottom=100*(a+b+c+d),width = .68, lw=0,color = color,alpha=.3,label = labels[4])
        g5 = plt.bar(position,100*f,bottom=100*(a+b+c+d+e),width = .68, lw=0,color = color,alpha=.1,label = labels[5])
        position+=1
        
    plt.ylabel('Genomic interactions (%)',fontsize=16)
    plt.ylim([0,100])
    plt.xlim([0.3,6.7])
    plt.title('{}\n'.format(hm_mark),fontsize=16,va='center')
    plt.legend([g5,g4,g3,g2,g1,g0],labels[::-1],loc=1,bbox_to_anchor=(1.32,1.05),fontsize=14,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,frameon=False)
    plt.axes().set_xticks([1.2,2.2,3.2,4.2,5.2,6.2]) 
    sns.despine(offset=0, trim=True)
    plt.axes().spines['bottom'].set_visible(False)
    plt.axes().set_xticklabels(ticklabels,rotation=45,ha='right',fontsize=15) 
    plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
#     plt.show()
    plt.savefig('{}/{}_genomic_distribution.pdf'.format(outdir,hm_mark),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)



# ==== main()

indir = 'f2_loop_genomic_feature'
outdir = 'f3_loop_genomic_distribution'
os.makedirs(outdir,exist_ok=True)

sum_df = pd.read_csv(outdir+os.sep+'loop_distribution_summary.csv',index_col=0)

# plot the figs
ticklabels = ['Vector','WT','$\Delta$cIDR','UTX_eIF$_{IDR}$','$\Delta$TPR']
hms = ['K27ACVEC','K27ACWT','K27ACDEL','K27ACEIF','K27ACTPR',]    
color = 'navy'
bar_plot(hms,ticklabels,'H3K27ac',outdir,color)   

hms = ['K4M3PCDH','K4M3WT','K4M3DEL3','K4M3EIF','K4M3DTPR',]
bar_plot(hms,ticklabels,'H3K4me3',outdir,color)   
    
