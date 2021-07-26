import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy import stats



def window_smooth(df,half_window):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df



def plot_composite_figs(df,outdir,basename,binding_type,norm_df):
    
    x = df.columns
    fig = plt.figure(figsize=(2.6,2.6))
    norm_factor=norm_df.loc[basename,'QC_dm6_f0x2']
    if binding_type=='rawCount':
        y = 1000000*df.mean().values/norm_factor
    else:
        y = df.mean().values
    plt.plot(x,y,color='red')
    
    # ==== ttest
    # s,p = stats.ttest_ind(y1.sum(axis=1).dropna().values,y4.sum(axis=1).dropna().values)
        
    plt.axes().set_xticks([x[0],x[int(len(x)/2)],x[-1]])
    plt.axvline(x=x[int(len(x)/2)],ls='--',lw=1,c='grey',zorder=0)
    plt.axes().set_xticklabels(['-2kb',0,'2kb'],rotation=0)
    plt.title(basename,fontsize=14)
    plt.ylabel('PROseq {}'.format(binding_type),fontsize=16)
    plt.xlabel('position relative to TSS',fontsize=16 )
    # plt.legend(fontsize=13,bbox_to_anchor=[1.72,1.05],loc="upper right",borderaxespad=0.5,labelspacing=.5,handletextpad=0.5,handlelength=1,frameon=False)
    plt.savefig(outdir+os.sep+'{}_{}.pdf'.format(binding_type,basename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    


        
        
## == main 

indir = 'f1_TSS_pattern'
outdir = 'f2_TSS_pattern_figs'
os.makedirs(outdir,exist_ok=True)

norm_df=pd.read_excel('../data_modules_revised/Proseq_norm_factor.xlsx',sheet_name='Normalization',index_col=0)
basenames= ['PCDH1PRO','PCDH2PRO','DEL31PRO','DEL32PRO','FL1PRO','FL2PRO',]
for basename in basenames[:2]:
    for binding_type in ['rawCount','RPKM',][:1]:
        csv_file='{}/{}_TSS_es2kb_bin200_{}.csv'.format(indir,basename,binding_type)    

        with open(csv_file) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=0)
        # print(df)
        # df = df[df.columns[75:126]]
        plot_composite_figs(df,outdir,basename,binding_type,norm_df)






