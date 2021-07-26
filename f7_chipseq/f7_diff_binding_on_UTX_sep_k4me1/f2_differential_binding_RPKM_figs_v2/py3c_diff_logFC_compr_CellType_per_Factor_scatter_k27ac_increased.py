import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde
from scipy import stats


cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}



compr_pairs = [['WT','Vector'],['DEL','WT'],['MT2','DEL']]
antibody_pairs = [['H3K27ac','H3K4me1'],['H3K27ac','H3K4me2'],\
                  ['H3K27ac','MLL4'],['UTX','H3K27ac'],\
                  ['UTX','H3K4me1'],['UTX','H3K4me2'],\
                  ['UTX','H3K4me3'],['UTX','MLL4'],\
                  ['H3K4me2','H3K4me1'],['H3K4me2','H3K4me3']]

factors = ['UTX','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
# peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
peak_files = ['UTX_peaks']

indir='../f0_data_integration_v3_RPKM/f3_combined_data_RPKM/'
outdir='f3c_diff_logFC_compr_CellType_per_Factor_scatter_k27ac_increased'
os.makedirs(outdir,exist_ok=True)

hm_log2fc_col='H3K27ac_WT_over_H3K27ac_Vector_log2FC'
hm_log2avg_col = 'H3K27ac_WT_over_H3K27ac_Vector_log2AVG'
fc_thre = 1.5
log2avg_thre = -2



stats_df = pd.DataFrame()
for peak_file in peak_files[:]:
    df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
    # for each factor, compr diff binding between cell types
    df = df[(df[hm_log2fc_col]> np.log2(fc_thre)) & (df[hm_log2avg_col]>log2avg_thre)] 
        
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
            
            if factor_logFC_y in df.columns:
                x = df[factor_logFC_x]
                y = df[factor_logFC_y]
                
                # == stats test
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
                s,p = stats.pearsonr(x, y)    
                output_prename = '{}_{}_binding_{}_over_{}_vs_{}_over_{}'.format(peak_file,factor,treatment_x,control_x,treatment_y,control_y)
                stats_df.loc[output_prename,'pearsonr_s'] = s
                stats_df.loc[output_prename,'pearsonr_p'] = p
                stats_df.loc[output_prename,'r_value'] = r_value
                stats_df.loc[output_prename,'p_value'] = p_value
    
                plt.figure(figsize=(2.3,2.3))
                data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
                z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
                z[np.where(np.isnan(z))] = 0.0
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
                # plt.scatter(x,y,rasterized=True,s=7,c='k')
                plt.axhline(y=0,c='gray',lw=.8,ls='--')
                plt.axvline(x=0,c='gray',lw=.8,ls='--')
                plt.xlabel('$\Delta${} ({} over {})'.format(factor[:7], cellType_labels[treatment_x],cellType_labels[control_x]),fontsize=13)
                plt.ylabel('$\Delta${} ({} over {})'.format(factor[:7], cellType_labels[treatment_y],cellType_labels[control_y]),fontsize=13)
                # plt.title('$\Delta${}'.format(factor))
                x_sort = np.sort(x)
                plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.9)
                # plt.text(.13,1.04,'$r={:.2f}, p={:.1e}$'.format(r_value,p_value),fontsize=13,transform=plt.axes().transAxes)
                plt.text(.97,.97,'$r={:.2f}$ '.format(r_value),fontsize=10,transform=plt.axes().transAxes,ha='right',va='top')
                plt.savefig(outdir+os.sep+'{}.pdf'.format(output_prename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
                plt.show()
                plt.close()
            
stats_df.to_csv(outdir+os.sep+'stats_summary.csv')   





