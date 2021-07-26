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
outdir='f2_diff_logFC_per_CellType_compr_Factor_scatter'
os.makedirs(outdir,exist_ok=True)


for peak_file in peak_files[:]:
    df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
    # for each pair of cell types, compr diff binding between two factors
    for compr_pair in compr_pairs[:1]:
        treatment = compr_pair[0]
        control = compr_pair[1]
        for antibody_pair in antibody_pairs[:]:
            factor1_logFC = '{}_{}_over_{}_{}_log2FC'.format(antibody_pair[0],treatment,antibody_pair[0],control)
            factor2_logFC = '{}_{}_over_{}_{}_log2FC'.format(antibody_pair[1],treatment,antibody_pair[1],control)
            x = df[factor1_logFC]
            y = df[factor2_logFC]

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
            plt.xlabel('$\Delta${}'.format(antibody_pair[0]))
            plt.ylabel('$\Delta${}'.format(antibody_pair[1]))
            plt.title('{} over {}'.format(compr_pair[0],compr_pair[1]))
            plt.savefig(outdir+os.sep+'{}_{}_vs_{}_compr_{}_over_{}.pdf'.format(peak_file,antibody_pair[0],antibody_pair[1],compr_pair[0],compr_pair[1]),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
            
    





