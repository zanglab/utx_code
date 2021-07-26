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




def return_deg(deseq_out,adjp=0.05,logfc=np.log2(1.5)):
    # read deg file, up genes
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_upgenes = deseq_out_upgenes.index
    # down genes
    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_dngenes = deseq_out_dngenes.index
    return deseq_out_upgenes,deseq_out_dngenes



def subplot_scatter(deg_df,dci_df,genes,color,rasterized=False):
    log2fc = deg_df.loc[genes,'log2FoldChange']
    dci_score = dci_df.loc[genes,'info']
    plt.scatter(log2fc,dci_score,color=color,s=5,rasterized=rasterized)



def plot_figs(deg_file,dci_file,outdir,hm,compr,compr_figname,dis_para):
    
    # read info
    deg_df = pd.read_csv(deg_file,index_col=0)
    dci_df = pd.read_csv(dci_file,sep='\t',index_col=4)
    # keep only those shared genes
    shared_genes= deg_df.index.intersection(dci_df.index)
    deg_df = deg_df.loc[shared_genes]
    dci_df = dci_df.loc[shared_genes]
    # select DEG
    upgenes,dngenes = return_deg(deg_df)
    
    # plot the fig
    plt.figure(figsize=(3,3))
    subplot_scatter(deg_df,dci_df,shared_genes,'grey',rasterized=True)
    subplot_scatter(deg_df,dci_df,upgenes,'r')
    subplot_scatter(deg_df,dci_df,dngenes,'b')
    
    # box_vals = []
    # for genes in [shared_genes,upgenes,dngenes]:
    #    box_vals.append(dci_df.loc[genes,'info'].values)
    # plt.violinplot(box_vals)
    
    plt.axhline(y=0,c='k',lw=1)
    plt.axvline(x=0,c='k',lw=1)
#     plt.xlim([-2,2])
    plt.title('{}'.format(compr))
    plt.xlabel('log2FoldChange')
    plt.ylabel('{} HiChIP DCI'.format(hm))
    plt.savefig(outdir+os.sep+'{}_{}_{}.pdf'.format(dis_para,hm,compr_figname),bbox_inches='tight',pad_inches=0.1,dpi=600)
#     plt.show()
    plt.close()
    
    pd.concat([deg_df,dci_df],axis=1).to_csv(outdir+os.sep+'_{}_{}_{}.csv'.format(dis_para,hm,compr_figname))
    
    

# ==== main() 

outdir = 'f1_expr_DCI_scatter'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)


for dis_para in ['dis200k','dis500k']:
    dci_dir='{}/f5_hichip/f1_hichip_bart3d/f1_promoter_DCI/promoter_DCI_{}'.format(project_dir,dis_para)

    # plot figs comparing each group of data
    compr='WT over Vector'
    compr_figname='WT_over_Vector'
    hm='H3K4me3'
    deg_file=expr_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv'
    dci_file=dci_dir+os.sep+'K4M3WT_over_K4M3PCDH_DCI_on_promoter.bed'
    plot_figs(deg_file,dci_file,outdir,hm,compr,compr_figname,dis_para)
    
    
    # 
    compr='WT over Vector'
    compr_figname='WT_over_Vector'
    hm='H3K27ac'
    deg_file=expr_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv'
    dci_file=dci_dir+os.sep+'K27ACWT_over_K27ACVEC_DCI_on_promoter.bed'
    plot_figs(deg_file,dci_file,outdir,hm,compr,compr_figname,dis_para)
    
    # 
    compr='$\Delta$cIDR over WT'
    compr_figname='DEL_over_WT'
    hm='H3K4me3'
    deg_file=expr_dir+os.sep+'treated_del_cIDR_vs_ctrl_WT.deseq2.csv'
    dci_file=dci_dir+os.sep+'K4M3DEL3_over_K4M3WT_DCI_on_promoter.bed'
    plot_figs(deg_file,dci_file,outdir,hm,compr,compr_figname,dis_para)
    
    # 
    compr='$\Delta$cIDR over WT'
    compr_figname='DEL_over_WT'
    hm='H3K27ac'
    deg_file=expr_dir+os.sep+'treated_del_cIDR_vs_ctrl_WT.deseq2.csv'
    dci_file=dci_dir+os.sep+'K27ACDEL_over_K27ACWT_DCI_on_promoter.bed'
    plot_figs(deg_file,dci_file,outdir,hm,compr,compr_figname,dis_para)
    
    



