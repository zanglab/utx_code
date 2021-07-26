import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")



def linear_regression(x,y):
    xmean = np.mean(x)
    ymean = np.mean(y)
    assert len(x)==len(y)
    sum1,sum2=0,0
    x_values,y_values = x.values,y.values
    for i in np.arange(len(x)):
        sum1 += (x_values[i]-xmean)*(y_values[i]-ymean)
        sum2 += np.power((x_values[i]-xmean),2)
    a = sum1/sum2
    b = ymean-a*xmean
    return a,b,xmean,ymean


    
    
def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ref(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*1.01 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),1)*1.01
    if compr_pos[1] ==3:
        y2 = y2-2
    if compr_pos[1] ==2:
        y2 = y2-2
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    # p_label='*'
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05:
        # if p<0.001:
            # p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=14)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=14)


def compr_plot(high_var_df,chip_marker,rpkm_pattern,outdir,):

    columns = ['WT', 'Vector', 'del_cIDR', 'del_TPR', 'UTX_eIFIDR',]
    colors = ['grey']*5
    positions = np.array([0,1,2,3,4])
    box_vals = []
    for i in positions:
        box_vals.append(np.sqrt(high_var_df[columns[i]].values))

    ## box with scatter
    plt.figure(figsize=(5,4))
    g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
                
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99)

    # for compr_pos in [[0,1,'b'],[0,2,'b'],[0,3,'b'],[0,4,'b']]:
    #     mark_pvalue(compr_pos,positions,box_vals)

    # plt.xlim([0,i+1])
    # plt.ylim([-.1,7])
    # plt.axes().set_xticks(positions+0.7)
    plt.axes().set_xticklabels(columns,rotation=30, ha='right',fontsize=17,color='k')
    plt.ylabel('sqrt(RPM) of {} at {} genes'.format(chip_marker,high_var_df.shape[0]),fontsize=17)
    plt.title(rpkm_pattern,fontsize=17)
    plt.show()
    figname=outdir+os.sep+'{}_{}.png'.format(chip_marker,rpkm_pattern)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()







# ==== main 
    
outdir = 'f4_compr_chip_binding_box'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
deg_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f2_deg/'.format(project_dir)
wt_up_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'
wt_dn_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'
up_genes = [i.strip() for i in open(wt_up_gene_file).readlines()]
dn_genes = [i.strip() for i in open(wt_dn_gene_file).readlines()]


chip_markers= ['H3K27ac','HA','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody','Macs2peak','Sicer2peak']
# celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody',]





for chip_marker in chip_markers[:]:
    for rpkm_pattern in rpkm_patterns[:]:
        csv_file='f1_highly_variable_chip_binding/{}_{}_RPM_residual.csv'.format(chip_marker,rpkm_pattern)
        # csv_file='../archived/f1_chip_binding_count_RPKM/f3_highly_variable_chip_binding//{}_{}_TPM_residual.csv'.format(chip_marker,rpkm_pattern)
        
        df = pd.read_csv(csv_file,index_col=0)
        # == select high-var genes
        diff_gene_df = df.loc[df.index.intersection(up_genes)]
        high_var_df = df[df['residual_pvalue']<0.05]
        compr_plot(diff_gene_df,chip_marker,rpkm_pattern,outdir,)
        
        
