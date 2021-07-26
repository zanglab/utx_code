import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from statistics import mode
    


# plot the figs

def loop_length_distribution(hm,cellTypes,indir,outdir,subdir):
    
    os.makedirs(outdir+os.sep+subdir,exist_ok=True)
    
    plt.figure(figsize=(2.6,2.6))
    all_length = []
    for celltype in cellTypes[:]:
        df = pd.read_csv('{}/{}/{}_{}.csv'.format(indir,subdir,hm,celltype),index_col=0)
        loop_length = [i.split('_') for i in df.index[:]]
        loop_length = np.abs([int(i[3])-int(i[1]) for i in loop_length])
        all_length = np.append(all_length,loop_length)
        kwargs = {'cumulative': False}
        g = sns.distplot(loop_length/1000,kde=True,norm_hist=False,hist=False,\
                          hist_kws=kwargs, kde_kws=kwargs, \
                      color = cellType_colors[celltype], label = cellType_labels[celltype],)
        x = g.lines[0].get_xdata() # Get the x data of the distribution
        y = g.lines[0].get_ydata() # Get the y data of the distribution
        mode_idx = np.argmax(y) # The id of the peak (maximum of y data)
        model_x,model_y = x[mode_idx],y[mode_idx]
        print(model_x,model_y)
        # plt.plot(model_x,model_y, '*', ms=13)

    plt.axvline(x=model_x,c='gray',lw=.8,ls='--')
    plt.text(model_x,0,' ${}$\n'.format(model_x.round(1)),fontsize=11,ha='left',va='center')
    plt.ylabel('Probability Density',fontsize=14)
    plt.text(.99,-.12,'kb',transform=plt.axes().transAxes)
    plt.xlim(xmin=0,xmax=500)
    plt.title('{}'.format(hm),fontsize=14)
    plt.legend(fontsize=12,borderaxespad=0.,labelspacing=.2,handletextpad=0.2,handlelength=1,frameon=False)
    plt.show()
    plt.savefig('{}/{}/{}.pdf'.format(outdir,subdir,hm),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    return model_x,model_y



# ==== main()
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


indir = '../f2_loop_anchor_annotation/f2_loop_genomic_feature/'
outdir = 'f1_loop_length_distribution_pdf'
subdirs=['data_1st_submission_rep_combined','data_1st_submission_sep_rep','data_202008']


# plot the figs
subdir=subdirs[0]
hm='H3K4me3'
cellTypes = ['Vector','WT','DEL','EIF']    
model_x,model_y = loop_length_distribution(hm,cellTypes,indir,outdir,subdir)   


# ==== batch2 H3K4me3
subdir=subdirs[2]
hm='H3K4me3'
cellTypes = ['Vector','WT','DEL','EIF','TPR']    
model_x,model_y = loop_length_distribution(hm,cellTypes,indir,outdir,subdir)   

# ==== batch2 H3K27ac
hm='H3K27ac'
model_x,model_y = loop_length_distribution(hm,cellTypes,indir,outdir,subdir)   




