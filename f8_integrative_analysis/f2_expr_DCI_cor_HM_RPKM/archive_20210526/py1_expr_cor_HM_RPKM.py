import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=11
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde



def plot_expression_cor_RPKM(df,stats_df,outdir,factor,count_type,compr_type):
    
    treatment,control = compr_type.split('_over_')
    expr_col = '{}_log2FoldChange'.format(compr_type)
    rpkm_col = '{}_{}_{}_log2FC'.format(count_type,factor,compr_type)
    x = df[expr_col]
    y = df[rpkm_col]
    
    # == stats test
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
    output_prename = '{}_{}_cor_expr_{}'.format(factor,count_type,compr_type)
    stats_df.loc[output_prename,'pearsonr_s'] = r_value
    stats_df.loc[output_prename,'pearsonr_p'] = p_value


    plt.figure(figsize=(2,2))
    data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
    z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
    z[np.where(np.isnan(z))] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
    # plt.scatter(x,y,rasterized=True,s=7,c='k')
    plt.axhline(y=0,c='gray',lw=.8,ls='--')
    plt.axvline(x=0,c='gray',lw=.8,ls='--')
    plt.xlabel('Expression log$_2$(fold change) \n {} over {}'.format(cellType_labels[treatment],cellType_labels[control]),fontsize=12)
    plt.ylabel('{} $\Delta${} \n {} over {}'.format(factor,count_type, cellType_labels[treatment],cellType_labels[control]),fontsize=12)
    # plt.title('$\Delta${}'.format(factor))
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.9)
    # plt.text(.13,1.04,'$r={:.2f}, p={:.1e}$'.format(r_value,p_value),fontsize=13,transform=plt.axes().transAxes)
    plt.text(.97,.97,'$r={:.2f}$'.format(r_value),fontsize=11,transform=plt.axes().transAxes,ha='right',va='top')
    plt.savefig(outdir+os.sep+'{}.pdf'.format(output_prename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()





cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

 
## ==== main
outdir = 'f1_expr_cor_HM_RPKM'
os.makedirs(outdir,exist_ok=True)

project_dir = "/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir = "/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
factors = ['H3K4me1','H3K4me3','H3K27ac','H3K27me3','H3K4me2APR']


combined_df = pd.read_csv('../data/gene_Expression_DCI_RP/fz_data_combined/TPM_DEseq2_DCI_RPKM_RP.csv',index_col=0)
# combined_df = combined_df.dropna()

# ==== plot the correlation between expression and RPKM/RP changes 

count_types = ['RP','RPKM']
compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']

stats_df = pd.DataFrame()
for factor in factors[:]:
    for count_type in count_types[:1]:
        for compr_type in compr_types[:]:
            plot_expression_cor_RPKM(combined_df,stats_df,outdir,factor,count_type,compr_type)

stats_df.to_csv(outdir+os.sep+'stats_summary.csv')   



