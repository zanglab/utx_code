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


    
def compr_plot(df1,df2,compr_col,marker_pair,rpkm_pattern,count_pattern,outdir,genes,color,label):
    
    plt.figure(figsize=(3,3))
    x,y = df1[compr_col],df2[compr_col]
    # x,y = -1*np.log10(df1[compr_col]),-1*np.log10(df2[compr_col])
    t = plt.scatter(x,y,s=3,c='darkgray',rasterized=True,alpha=1,label='others')

    genes = df1.index.intersection(genes)
#     dn_genes = df1.index.intersection(dn_genes)
    t = plt.scatter(x[genes],y[genes],s=4,c=color,rasterized=True,alpha=1,label=label)
#     t = plt.scatter(x[dn_genes],y[dn_genes],s=4,c='r',rasterized=True,alpha=1,label='Vector up')


    # regression
    # slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    # x_sort = np.sort(x)
    # plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    # plt.text(.55,1.03,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=15,transform=plt.axes().transAxes)
    plt.legend(markerscale=3,fontsize=10,frameon=False)
    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.title('{}'.format(rpkm_pattern),fontsize=16)
    plt.xlabel('{} residual'.format(marker_pair[0]),fontsize=16)
    plt.ylabel('{} residual'.format(marker_pair[1]),fontsize=16)
    figname = '{}/{}_vs_{}_{}_{}_{}.png'.format(outdir,marker_pair[0],marker_pair[1],rpkm_pattern,count_pattern,'_'.join(label.split(' ')))
#     plt.show()
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()

    
    


# ==== main 
    
outdir = 'f2_compr_chip_binding_residual_figs'
os.makedirs(outdir,exist_ok=True)
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
deg_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f2_deg/'.format(project_dir)
wt_up_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_upgenes.txt'
wt_dn_gene_file = deg_dir+os.sep+'treated_WT_vs_ctrl_Vector.deseq2.csv_adjp0.05_logfc0.32_dngenes.txt'
up_genes = [i.strip() for i in open(wt_up_gene_file).readlines()]
dn_genes = [i.strip() for i in open(wt_dn_gene_file).readlines()]


chip_markers= ['HA','H3K27ac','H3K27me3','H3K4me1_rep','H3K4me1','H3K4me3','MLL4','MLL4SC','UTX']
# rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody','Macs2peak','Sicer2peak']
# celltypes = ['del_cIDR','del_TPR','UTX_eIFIDR','WT','Vector']

marker_pairs=[['H3K27ac','H3K4me3'],['H3K27ac','H3K4me1'],['H3K27ac','H3K4me1_rep'],['H3K27ac','H3K27me3'],['H3K27ac','MLL4SC'],]
rpkm_patterns = ['Promoter','Promoter_expand_10kb','GeneBody']
count_patterns = ['RPKM','RPM']

for marker_pair in marker_pairs[:]:
    for rpkm_pattern in rpkm_patterns[:]:
        for count_pattern in count_patterns:
            csv_file1='f1_highly_variable_chip_binding/{}_{}_{}_residual.csv'.format(marker_pair[0],rpkm_pattern,count_pattern)
            csv_file2='f1_highly_variable_chip_binding/{}_{}_{}_residual.csv'.format(marker_pair[1],rpkm_pattern,count_pattern)
            df1 = pd.read_csv(csv_file1,index_col=0)
            df2 = pd.read_csv(csv_file2,index_col=0)
            shared_index = df1.index.intersection(df2.index)
            df1 = df1.loc[shared_index]
            df2 = df2.loc[shared_index]
        
            # compr_col='residual_pvalue'
            compr_col='residual'
            compr_plot(df1,df2,compr_col,marker_pair,rpkm_pattern,count_pattern,outdir,up_genes,'royalblue','Vector down')
            compr_plot(df1,df2,compr_col,marker_pair,rpkm_pattern,count_pattern,outdir,dn_genes,'r','Vector up')
        
        
            fisher_p_thre = 0.05
            sig1=df1[df1['residual_pvalue']<fisher_p_thre]
            sig2=df2[df2['residual_pvalue']<fisher_p_thre]
            shared_sig_genes = sig1.index.intersection(sig2.index)
            shared_sig_genes_down = shared_sig_genes.intersection(up_genes)
            a,b = len(shared_sig_genes_down),len(shared_sig_genes)-len(shared_sig_genes_down)
            c,d = len(shared_sig_genes),len(shared_index)-len(shared_sig_genes)
            s,p = stats.fisher_exact([[a,b],[c,d]])
            print(marker_pair,rpkm_pattern,count_pattern,a,b,c,d, s,p)
       
