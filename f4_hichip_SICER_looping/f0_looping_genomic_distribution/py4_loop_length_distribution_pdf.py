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
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

    


# plot the figs

def loop_length_distribution(hms,cellTypes,colors,indir,outdir,flag):
    
    plt.figure(figsize=(3,2.6))
    
    for ii in np.arange(len(hms)):
        hm, cellType, color = hms[ii], cellTypes[ii], colors[ii]
        df = pd.read_csv('{}/{}.csv'.format(indir,hm),index_col=0)
        loop_length = df['start2']-df['start1']
        ax = sns.distplot(loop_length/1000,kde=True,norm_hist=False,hist=False,\
                      color = color, label = cellType,)
        
    plt.ylabel('PDF',fontsize=16)
    plt.text(1.05,-.125,'kb',transform=plt.axes().transAxes)
    plt.xlim(xmin=-12,xmax=500)
    plt.title('{}'.format(flag),fontsize=16)
    plt.legend(fontsize=13,borderaxespad=0.,labelspacing=.2,handletextpad=0.2,handlelength=1,frameon=False)
    plt.show()
    plt.savefig('{}/{}.pdf'.format(outdir,flag),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)



# ==== main()

indir = 'f2_loop_genomic_feature'
outdir = 'f4_loop_length_distribution_pdf'
os.makedirs(outdir,exist_ok=True)

# list
cellTypes = ['Vector','WT','$\Delta$cIDR','UTX_eIF$_{IDR}$','$\Delta$TPR']
colors = ['tab:blue','tab:red','k','tab:purple','tab:green',]

hms = ['K27ACVEC','K27ACWT','K27ACDEL','K27ACEIF','K27ACTPR',]    
loop_length_distribution(hms,cellTypes,colors,indir,outdir,'H3K27ac',)   

hms = ['K4M3PCDH','K4M3WT','K4M3DEL3','K4M3EIF','K4M3DTPR',]
loop_length_distribution(hms,cellTypes,colors,indir,outdir,'H3K4me3',)   

    
