import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn2_plot(deg_list_t,deg_list_c,la,lb,figname,color1='r',color2='k'):

    a = set([i.strip() for i in open(deg_list_t).readlines()])
    b = set([i.strip() for i in open(deg_list_c).readlines()])

    # fisher exact test
    total=20794
    va = len(a.intersection(b))
    vb = len(a) - len(a.intersection(b))
    vc = len(b) - len(a.intersection(b))
    vd = total - len(a) - len(b) + len(a.intersection(b))
    s,p = stats.fisher_exact([[va,vb],[vc,vd]])
    print(len(a),len(b))
    print('odds ratio={:.2f}, pvalue={:.2e}'.format(s,p))
    
    # plot venn diagram
    plt.figure(figsize=(4,4))
    out = venn2([a,b],set_labels=(la,lb),set_colors=(color1,color2), alpha=0.5)
    
    for text in out.set_labels:
        text.set_fontsize(16)
    for text in out.subset_labels:
        try:
            text.set_fontsize(16)
        except:
            pass

    if s>0 and p<0.05:
        plt.text(x=.2,y=-.25,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches=0.1,transparent=True)
#     plt.show()
    plt.close()
    
    with open(figname+'.txt','w') as outf:
        outf.write('\n'.join(a.intersection(b))+'\n')




if 1:
    
    outdir = 'f7_deg_venn'
    os.makedirs(outdir,exist_ok=True)

    file_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
#     file_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
    indir = '{}/f6_deg/f2_deg'.format(file_dir)
    
    # WT up DEL up
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    deg_list_c = '{}/treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-up_DEL-up_venn.png'
    la,lb='up-genes in \nWT over Vector','up-genes in \nDEL over WT'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname)
    
    # WT up DEL down
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    deg_list_c = '{}/treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-up_DEL-down_venn.png'
    la,lb='up-genes in \nWT over Vector','down-genes in \nDEL over WT'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname)
    
    # WT down DEL up
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    deg_list_c = '{}/treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-dn_DEL-up_venn.png'
    la,lb='down-genes in \nWT over Vector','up-genes in \nDEL over WT'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname)
    
    # WT down DEL down
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    deg_list_c = '{}/treated_del_cIDR_vs_ctrl_WT.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-dn_DEL-down_venn.png'
    la,lb='down-genes in \nWT over Vector','down-genes in \nDEL over WT'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname)


    # WT up EIF up
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    deg_list_c = '{}/treated_UTX_eIFIDR_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-up_EIF-up_venn.png'
    la,lb='up-genes in \nWT over Vector','up-genes in \nEIF over Vector'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname,color2='goldenrod')

    # WT up EIF down
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    deg_list_c = '{}/treated_UTX_eIFIDR_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-up_EIF-down_venn.png'
    la,lb='up-genes in \nWT over Vector','down-genes in \nEIF over Vector'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname,color2='goldenrod')
    
    # WT down EIF up
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    deg_list_c = '{}/treated_UTX_eIFIDR_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-dn_EIF-up_venn.png'
    la,lb='down-genes in \nWT over Vector','up-genes in \nEIF over Vector'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname,color2='goldenrod')
   
    # WT down EIF down
    deg_list_t = '{}/treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    deg_list_c = '{}/treated_UTX_eIFIDR_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'.format(indir)
    figname = outdir+os.sep+'WT-dn_EIF-down_venn.png'
    la,lb='down-genes in \nWT over Vector','down-genes in \nEIF over Vector'
    venn2_plot(deg_list_t,deg_list_c,la,lb,figname,color2='goldenrod')

