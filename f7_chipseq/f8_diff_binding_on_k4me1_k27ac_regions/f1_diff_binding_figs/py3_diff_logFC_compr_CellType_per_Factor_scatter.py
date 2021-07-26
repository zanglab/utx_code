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
from scipy import stats



compr_pairs = [['WT','Vector'],['DEL','WT'],['EIF','DEL']]
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
peak_files = ['H3K27ac_peaks','H3K27ac_islands_increased',
              'H3K4me1_peaks','H3K4me1_islands_increased',
              'MLL4_peaks','MLL4_islands_increased']

indir='../f0_data_integration/f2_combined_data/'
outdir='f3_diff_logFC_compr_CellType_per_Factor_scatter'
os.makedirs(outdir,exist_ok=True)

stats_df = pd.DataFrame()
for peak_file in peak_files[:]:
    df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
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

            plt.figure(figsize=(3,3))
            data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
            z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
            z[np.where(np.isnan(z))] = 0.0
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
            # plt.scatter(x,y,rasterized=True,s=7,c='k')
            plt.axhline(y=0,c='gray',lw=.8,ls='--')
            plt.axvline(x=0,c='gray',lw=.8,ls='--')
            plt.xlabel('$\Delta${} ({} over {})'.format(factor, treatment_x,control_x))
            plt.ylabel('$\Delta${} ({} over {})'.format(factor, treatment_y,control_y))
            # plt.title('$\Delta${}'.format(factor))
            x_sort = np.sort(x)
            plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.9)
            plt.text(.13,1.04,'$r={:.2f}, p={:.1e}$'.format(r_value,p_value),fontsize=13,transform=plt.axes().transAxes)
            # plt.text(.50,1.03,'$P={:.1e}$ '.format(p_value),fontsize=13,transform=plt.axes().transAxes)
            plt.savefig(outdir+os.sep+'{}.pdf'.format(output_prename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
            
stats_df.to_csv(outdir+os.sep+'stats_summary.csv')   





