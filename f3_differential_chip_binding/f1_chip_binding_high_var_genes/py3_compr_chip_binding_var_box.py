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
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")




def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),96)*1.01 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    y = y + compr_pos[1]*.7
    p_label='*'
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    
    if s<0:
        p_label = p_label+'(+)'
    else:
        p_label = p_label+'(-)'
        
    if p<0.05:
        # if p<0.001:
            # p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color='k',fontsize=13)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color='k',fontsize=13)


    
def box_plot(box_vals,positions,chip_marker,rpkm_pattern,figname):  

    plt.figure(figsize=(3,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    
    colors = ['grey','r','b']            
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99)

    positions=[0,1,2]
    for compr_pos in [[0,1,'t'],[0,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)
    
    # plt.xlim([0,i+1])
    # plt.ylim([-.1,7])
    # plt.axes().set_xticks(positions+0.7)
    plt.axes().set_xticklabels(['All genes','Vector Up','Vector Down'],rotation=30, ha='right',fontsize=16,color='k')
    # plt.legend([g["boxes"][0],g["boxes"][1],g["boxes"][2]],['Control','ENCODE','GTEx'],fontsize=16,borderaxespad=0.1,labelspacing=.1,handletextpad=0.2,loc="upper left",frameon=False)
    # plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')    
    plt.ylabel('{} residual \n at {}'.format(chip_marker,rpkm_pattern),fontsize=16)
#     plt.show()
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()




# ==== main 
    
outdir = 'f3_compr_chip_binding_residual_figs_box'
os.makedirs(outdir,exist_ok=True)
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
deg_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f2_deg/'.format(project_dir)
wt_up_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'
wt_dn_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'
up_genes = [i.strip() for i in open(wt_up_gene_file).readlines()]
dn_genes = [i.strip() for i in open(wt_dn_gene_file).readlines()]


chip_markers= ['H3K27ac','HA','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']
# rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody','Macs2peak','Sicer2peak']
# celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']

marker_pairs=[['H3K27ac','H3K4me3'],['H3K27ac','H3K4me1'],['H3K27ac','H3K4me1_rep'],['H3K27ac','H3K27me3'],['H3K27ac','MLL4SC'],]
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody']
count_patterns = ['RPKM','RPM']
compr_col = 'residual'

for chip_marker in chip_markers[:]:
# for marker_pair in marker_pairs[:1]:
    for rpkm_pattern in rpkm_patterns[:]:
        for count_pattern in count_patterns[:1]:
            csv_file='f1_highly_variable_chip_binding/{}_{}_{}_residual.csv'.format(chip_marker,rpkm_pattern,count_pattern)
            df = pd.read_csv(csv_file,index_col=0)
      
            # compr_col='residual_pvalue'
            # compr_col='residual'
            # compr_plot(df1,df2,compr_col,marker_pair,rpkm_pattern,count_pattern,outdir,up_genes,'royalblue','Vector down')
            # compr_plot(df1,df2,compr_col,marker_pair,rpkm_pattern,count_pattern,outdir,dn_genes,'r','Vector up')
        
            fisher_p_thre = 0.1
            # df = df[df['residual_pvalue']<fisher_p_thre]
            wt_up_df = df.loc[df.index.intersection(up_genes)]
            wt_dn_df = df.loc[df.index.intersection(dn_genes)]
            
            box_vals = [df[compr_col],wt_dn_df[compr_col],wt_up_df[compr_col]]
            positions = [0,1,2]
            figname = '{}/{}_{}_{}.png'.format(outdir,chip_marker,rpkm_pattern,count_pattern)
            box_plot(box_vals,positions,chip_marker,rpkm_pattern,figname)
                
                
            # box_vals
            # sig=df[df['residual_pvalue']<fisher_p_thre]
            # sig2=df2[df2['residual_pvalue']<fisher_p_thre]
            # shared_sig_genes = sig1.index.intersection(sig2.index)
            # shared_sig_genes_down = sig.index.intersection(up_genes)
            # a,b = len(shared_sig_genes_down),len(shared_sig_genes)-len(shared_sig_genes_down)
            # c,d = len(shared_sig_genes),len(shared_index)-len(shared_sig_genes)
            # s,p = stats.fisher_exact([[a,b],[c,d]])
            # print(marker_pair,rpkm_pattern,count_pattern, s,p)
       
