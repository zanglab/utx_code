import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
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





def return_differential_values(indir,antibody,genomic_region,compr_pair):
    csv_file='{}/{}_{}_differential_RPKM.csv'.format(indir,antibody,genomic_region)
    df = pd.read_csv(csv_file,index_col=0)
    diff_df = df[[compr_pair]]
    # thre = np.percentile(np.abs(diff_df),80)
    # diff_df = diff_df[np.abs(diff_df)>thre]
    return diff_df



project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# cellTypes = ['Vector','WT','DEL','EIF']
# antibodys = ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
# antibodys = ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']

genomic_regions = ['UDHS','Promoter']
compr_pairs = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']
antibody_pairs = [['H3K27ac','H3K4me1'],['H3K27ac','H3K4me2'],\
                  ['H3K27ac','MLL4'],['UTX','H3K27ac'],\
                  ['UTX','H3K4me1'],['UTX','H3K4me2'],\
                  ['UTX','H3K4me3'],['UTX','MLL4'],\
                  ['H3K4me2','H3K4me1'],['H3K4me2','H3K4me3']]

# indir='f1_differential_RPKM_MA'
indir='f1b_differential_RPKM_MA'
outdir='f2_differential_RPKM_compr_antibody_scatter'
os.makedirs(outdir,exist_ok=True)

# peak_overlapped_UDHS_file = '{}/f7_chipseq/data_peak_overlap/overlapped/udhs_overlap_ALL.bed'.format(project_dir)
# peak_overlapped_UDHS = pd.read_csv(peak_overlapped_UDHS_file,sep='\t',index_col=3,header=None)

# for celltype in cellTypes:
for antibody_pair in antibody_pairs[:]:
    for genomic_region in genomic_regions[:]:
        for compr_pair in compr_pairs[:]:
            diff_score1 = return_differential_values(indir,antibody_pair[0],genomic_region,compr_pair)
            diff_score2 = return_differential_values(indir,antibody_pair[1],genomic_region,compr_pair)
            # x = diff_score1.loc[peak_overlapped_UDHS.index]
            # y = diff_score2.loc[peak_overlapped_UDHS.index]
            shared_index = diff_score1.index.intersection(diff_score2.index)
            x = diff_score1.loc[shared_index][compr_pair].values
            y = diff_score2.loc[shared_index][compr_pair].values

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
            plt.xlabel('$\Delta${} at {}'.format(antibody_pair[0],genomic_region))
            plt.ylabel('$\Delta${} at {}'.format(antibody_pair[1],genomic_region))
            plt.title('{}'.format(compr_pair,))
            plt.savefig(outdir+os.sep+'{}_{}_vs_{}_compr_{}.pdf'.format(genomic_region,antibody_pair[0],antibody_pair[1],compr_pair),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
       






