import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn



# ==== main() 

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}


    
project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
master_file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/f0_data_integration/f2_combined_data'.format(project_dir)
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)

outdir = 'f3_utx_binding_DEG_logFC_scatter'
os.makedirs(outdir,exist_ok=True)

peak_files = ['UTX_peaks','UTX_islands','UTXFEB_islands','UTXFEB_peaks']
k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
fc_thres = [1.5,2]
fc_thre = 1.5
log2avg_thre = 0

extend_dis = [0,2000,10000,50000]
compr_types = [['WT_over_Vector','DEL_over_WT'],['DEL_over_WT','EIF_over_DEL']]

# expression of genes near UTX           
for peak_file in peak_files[:]:
    # read the master dataframe
    master_df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(master_file_dir,peak_file),index_col=0)
    master_df_tmp = master_df[(master_df[k4me1_log2fc_col]> np.log2(fc_thre)) & (master_df[k4me1_log2avg_col]>log2avg_thre)] 
    for compr_type in compr_types[:]:
        compr_x =  compr_type[0]
        compr_y = compr_type[1]
        # log2FC of all genes
        all_log2FC_x =  deg_df['{}_log2FoldChange'.format(compr_x)].dropna()
        all_log2FC_y =  deg_df['{}_log2FoldChange'.format(compr_y)].dropna()
        # UTX nearby genes
        for dis in extend_dis[:]:
            utx_file = 'f1_UTX_binding_promoter_overlap/{}_tss_es{}bp.csv'.format(peak_file,dis)
            utx_df = pd.read_csv(utx_file,sep='\t',index_col=4)
            utx_df = utx_df.loc[master_df.index].dropna()
            overlapped_genes = utx_df[utx_df['IfOverlap']==1]['overlapped_genes']
            overlapped_genes = [gene for ele in overlapped_genes for gene in ele.split(',')]
            overlapped_genes = set(overlapped_genes).intersection(deg_df.index)
            # logFC of overlapped genes
            overlapped_log2FC_x = deg_df.loc[overlapped_genes]['{}_log2FoldChange'.format(compr_x)].dropna()
            overlapped_log2FC_y = deg_df.loc[overlapped_genes]['{}_log2FoldChange'.format(compr_y)].dropna()
            
            # scatter plot
            plt.figure(figsize=(3,3))
            plt.scatter(all_log2FC_x,all_log2FC_y,c='grey',s=5,rasterized=True,label='All genes')

            plt.scatter(overlapped_log2FC_x,overlapped_log2FC_y,c='red',s=5,rasterized=True,label='w/ UTX')
#             x,y = overlapped_log2FC_x,overlapped_log2FC_y
#             data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
#             z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
#             z[np.where(np.isnan(z))] = 0.0
#             idx = z.argsort()
#             x, y, z = x[idx], y[idx], z[idx]
#             g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=5,rasterized=True,label='w/ UTX')
            plt.axhline(y=0,c='k',lw=1)
            plt.axvline(x=0,c='k',lw=1)
            # treatment,control = compr_type.split('_')[0], compr_type.split('_')[-1]
            # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
            plt.legend(fontsize=13,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,\
                       handlelength=1,loc="upper right",markerscale=3)
            plt.xlabel('log2FC of {}'.format(compr_x))
            plt.ylabel('log2FC of {}'.format(compr_y))
            plt.savefig('{}/{}_dis{}bp_{}_vs_{}.png'.format(outdir,peak_file,dis,compr_x,compr_y),bbox_inches='tight',pad_inches=0.1,dpi=600)
            plt.show()
            plt.close()

