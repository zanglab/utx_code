import pandas as pd
import numpy as np
import glob,os,sys
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]



def plot_figs(outdir,subdir,basename,inf):

    inf  = inf.iloc[:10,:]
    y = [-np.log10(i) for i in inf["Adjusted P-value"]]
    if len(y)==0:
        return

    fig, ax = plt.subplots(figsize=(2.8,len(inf.index)*0.35))  
    plt.barh(np.arange(len(inf)),y,align='center', alpha=0.9,height=0.7,color='salmon')
    plt.yticks(np.arange(len(inf)),[i for i in inf.index])
    #plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
    plt.xlim([0,max(y)*1.2])
    plt.ylim([-0.8,len(y)-0.2])
    plt.gca().invert_yaxis()
    #plt.title(infile,verticalalignment='bottom')
    plt.xlabel("-log$_{{10}}$ (adj. $P$)",fontsize=18,verticalalignment='top')
    for i,v in enumerate(y):
        ax.text(v + 0.03*y[0],i +0.4,str('{:.1f}'.format(v)),fontsize=13)
        

    plt.title('{} {}\n{}'.format(subdir.split('.')[0],subdir.split('_')[-1],basename),fontsize=15)
#     plt.title('{} {} sites'.format(cancertype.split('_')[0],bindingtype))
    plt.savefig(outdir+os.sep+'{}_enrichr.pdf'.format(basename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    




# ==== main

indir = 'f3_enrichr'
outdir = 'f4_enrichr_figs'
os.makedirs(outdir,exist_ok=True)

subdirs = ['treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_dngenes',
'treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_upgenes',
'treated_UTX_eIFIDR_vs_ctrl_del_cIDR.deseq2.csv_adjp0.05_logfc0.32_dngenes',
'treated_UTX_eIFIDR_vs_ctrl_del_cIDR.deseq2.csv_adjp0.05_logfc0.32_upgenes',
'treated_UTX_eIFIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_dngenes',
'treated_UTX_eIFIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_upgenes',
'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes',
'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes']


for subdir in subdirs:
    # sub dir for each pathway
    sub_outdir = outdir+os.sep+subdir
    os.makedirs(sub_outdir,exist_ok=True)
    # plot each fig
    result_files = glob.glob(indir+os.sep+subdir+os.sep+'*.txt')
    for result_file in result_files:
        basename = os.path.basename(result_file).split('.')[0]
        print(subdir,basename)
        with open(result_file) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=0)
        df = df.sort_values(by=['Adjusted P-value'],ascending=True)  
        df = df[df['Adjusted P-value']<0.05]
        plot_figs(sub_outdir,subdir,basename,df)
            
            
            
