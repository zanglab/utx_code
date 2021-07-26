import sys,argparse
import os,glob
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
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from scipy.interpolate import interpn
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def scatter_all_chip_binding(df,columns,figname,xlabel,ylabel):

    fig = plt.figure(figsize=(4,4))
    x,y = np.log2(df[columns[0]]+1),df[columns[1]]
    g = plt.scatter(x,y,c='lightgrey',s=1,marker='o')

    df1 = df[df['padj']<0.05]
    x,y = np.log2(df1[columns[0]]+1),df1[columns[1]]
    g = plt.scatter(x,y,c='red',s=1,marker='o')
    plt.axhline(y=0.42,lw=2,ls='--',c='k')
    plt.axhline(y=-0.42,lw=2,ls='--',c='k')
    plt.xlabel(xlabel,fontsize=15)
    plt.ylabel(ylabel,fontsize=15)
#     plt.ylim([-0.9,0.9])
    plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()


def scatter_all_chip_binding_2(df,columns,figname,xlabel,ylabel,adjp,logfc,basename):

    fig = plt.figure(figsize=(4,4))
    x,y = df[columns[0]],-1*np.log10(df[columns[1]])
    g = plt.scatter(x,y,c='grey',s=1,marker='o')

    # up-regulated genes
    df1 = df[(df['log2FoldChange']>logfc) &(df['padj']<adjp)]
    x,y = df1[columns[0]],-1*np.log10(df1[columns[1]])
    g = plt.scatter(x,y,c='red',s=1,marker='o')

    # down-regulated genes
    df1 = df[(df['log2FoldChange']<-1*logfc) &(df['padj']<adjp)]
    x,y = df1[columns[0]],-1*np.log10(df1[columns[1]])
    g = plt.scatter(x,y,c='blue',s=1,marker='o')

    plt.axvline(x=0,lw=1,ls='--',c='k')
    plt.axhline(y=-1*np.log10(adjp),lw=1,ls='--',c='k')
    plt.xlabel(xlabel,fontsize=15)
    plt.ylabel(ylabel,fontsize=15)
    plt.title(basename+'\n')
    # plt.ylim([-.1,20])
    plt.xlim([-2.5,2.5])
    plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
#     plt.show()
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()



def return_diff_genes_deseq2(csv_file,adjp,logfc):

    deseq_out = pd.read_csv(csv_file,index_col=0)
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    # deseq_out_upgenes = set(deseq_out_upgenes['GeneID'])

    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    # deseq_out_dngenes = set(deseq_out_dngenes['GeneID'])  
    return deseq_out_upgenes,deseq_out_dngenes



# == main ==
if 1:
    outdir = 'f6_deg_scatter'
    os.makedirs(outdir,exist_ok=True)

    file_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
    file_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
    deseq2_file1 = '{}/f6_deg/f1_deseq2_out/treated_WT_vs_ctrl_Vector.deseq2.csv'.format(file_dir)
    deseq2_file2 = '{}/f6_deg/f1_deseq2_out/treated_del_cIDR_vs_ctrl_WT.deseq2.csv'.format(file_dir)
    adjp,logfc = 0.05,np.log2(1.25)
    # upgenes,dngenes = return_diff_genes_deseq2(deseq2_file,adjp,logfc)
    # deg = pd.concat([upgenes,dngenes]).sort_values(by=['log2FoldChange'],ascending=False)
    
    df1 = pd.read_csv(deseq2_file1,sep=',',index_col=0)
    df1 = df1[(np.abs(df1['log2FoldChange'])>1*logfc) &(df1['padj']<adjp)]
    
    df2 = pd.read_csv(deseq2_file2,sep=',',index_col=0)
    df2 = df2.loc[df1.index]
    
    x,y = df1['log2FoldChange'],df2['log2FoldChange']
    # x,y = np.log2(df1['baseMean']),df1['log2FoldChange']

    fig = plt.figure(figsize=(4,4))
    g = plt.scatter(x,y,c='grey',s=1,marker='o')

    plt.axvline(x=0,lw=1,ls='--',c='k')
    plt.axhline(y=0,lw=1,ls='--',c='k')
    # plt.xlabel(xlabel,fontsize=15)
    # plt.ylabel(ylabel,fontsize=15)
    # plt.title(basename+'\n')
    # plt.ylim([-.1,20])
    # plt.xlim([-2.5,2.5])
    plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    plt.show()
    # plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()
