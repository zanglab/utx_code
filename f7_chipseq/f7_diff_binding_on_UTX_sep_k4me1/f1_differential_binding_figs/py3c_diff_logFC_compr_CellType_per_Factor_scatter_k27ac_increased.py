import sys,argparse
import os,glob
import numpy as np
import pandas as pd
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
from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde





compr_pairs = [['WT','Vector'],['DEL','WT'],['EIF','DEL']]
antibody_pairs = [['H3K27ac','H3K4me1'],['H3K27ac','H3K4me2'],\
                  ['H3K27ac','MLL4'],['UTX','H3K27ac'],\
                  ['UTX','H3K4me1'],['UTX','H3K4me2'],\
                  ['UTX','H3K4me3'],['UTX','MLL4'],\
                  ['H3K4me2','H3K4me1'],['H3K4me2','H3K4me3']]

factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']

indir='../f0_data_integration/f2_combined_data/'
outdir='f3c_diff_logFC_compr_CellType_per_Factor_scatter_k27ac_increased'
os.makedirs(outdir,exist_ok=True)

# k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
# k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
hm_log2fc_col='H3K27ac_WT_over_H3K27ac_Vector_log2FC'
hm_log2avg_col = 'H3K27ac_WT_over_H3K27ac_Vector_log2AVG'
fc_thres = [2,1.5]
log2avg_thre = 0

for peak_file in peak_files[:]:
    df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
    for fc_thre in fc_thres:
        df_tmp = df[(df[hm_log2fc_col]> np.log2(fc_thre)) & (df[hm_log2avg_col]>log2avg_thre)] 
        print(peak_file, fc_thre, df_tmp.shape)
        # for each factor, compr diff binding between cell types
        for factor in factors[:]:
            for ii in np.arange(len(compr_pairs)-1)[:]:
                # compr_pairs[ii]
                treatment_x = compr_pairs[ii][0]
                control_x = compr_pairs[ii][1]
                # compr_pairs[ii+1]
                treatment_y = compr_pairs[ii+1][0]
                control_y = compr_pairs[ii+1][1]                
                factor_logFC_x = '{}_{}_over_{}_{}_log2FC'.format(factor,treatment_x,factor,control_x)
                factor_logFC_y = '{}_{}_over_{}_{}_log2FC'.format(factor,treatment_y,factor,control_y)
            
                x = df_tmp[factor_logFC_x]
                y = df_tmp[factor_logFC_y]
    
                plt.figure(figsize=(3,3))
                data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
                z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
                z[np.where(np.isnan(z))] = 0.0
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
                # plt.scatter(x,y,rasterized=True,s=7,c='k')
                plt.axhline(y=0,c='gray',lw=1,ls='--')
                plt.axvline(x=0,c='gray',lw=1,ls='--')
                plt.xlabel('$\Delta${} ({} over {})'.format(factor, treatment_x,control_x))
                plt.ylabel('$\Delta${} ({} over {})'.format(factor, treatment_y,control_y))
                # plt.title('$\Delta${}'.format(factor))
                plt.savefig(outdir+os.sep+'fcthre_{}_{}_{}_binding_{}_over_{}_vs_{}_over_{}.pdf'.format(fc_thre,peak_file,factor,treatment_x,control_x,treatment_y,control_y),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
                plt.show()
                plt.close()
                
        





