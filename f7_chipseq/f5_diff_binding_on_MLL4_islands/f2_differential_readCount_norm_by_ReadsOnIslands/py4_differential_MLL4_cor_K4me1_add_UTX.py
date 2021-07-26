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
from scipy.stats import gaussian_kde





def return_differential_values(indir,antibody,compr_pair,norm_col,genomic_region):
#     if antibody=='UTXFEB' or antibody=='UTX':
#         genomic_region='_es500bp'
#         norm_col='total'
#     else:
#         genomic_region='_es2kb'
#         norm_col='total_in_islads'
        
    csv_file='{}/_{}_{}_{}_differential_RPKM.csv'.format(indir,norm_col,antibody,genomic_region)
    df = pd.read_csv(csv_file,index_col=0)
    # df = df[df[df.columns[:4]].mean(axis=1)>2]
    diff_df = df[[compr_pair]]
    # thre = np.percentile(np.abs(diff_df),80)
    # diff_df = diff_df[np.abs(diff_df)>thre]
    return diff_df,genomic_region



project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

genomic_regions = ['_es2kb','_es500bp','']
compr_pairs = ['WT_over_Vector_log2FC','DEL_over_WT_log2FC','EIF_over_DEL_log2FC']
antibody_pairs = [['MLL4','H3K4me1'],['MLL4','H3K4me2']]
norm_cols = ['total_in_islads','total']

MLL4_UTX_overlap_dir='{}/f7_chipseq/f5_diff_binding_on_MLL4_islands/data_UTX_overlapped_MLL4/overlapped'.format(project_dir)
overlap_types=['202011UTX_WT_over_Vector_islands','202102UTX_WT_over_Vector_islands',\
               '202011UTX_WT_over_Vector_peaks','202102UTX_WT_over_Vector_peaks']

indir='f1_differential_binding_MA'
outdir='f4_differential_MLL4_cor_K4me1'
os.makedirs(outdir,exist_ok=True)


# for celltype in cellTypes:
for overlap_type in overlap_types[:]:
    utx_overlapped_MLL4_file='{}/MLL4_on_{}.bed'.format(MLL4_UTX_overlap_dir,overlap_type)
    utx_overlapped_MLL4_df = pd.read_csv(utx_overlapped_MLL4_file,sep='\t',header=None)
    overlapped_eles=utx_overlapped_MLL4_df.loc[:,:2].astype(str).values
    overlapped_ids = ['_'.join(ii) for ii in overlapped_eles]
    # then check the diff MLL4/k4me1 binding, highlightin thouse UTX binding sites
    for antibody_pair in antibody_pairs[:]:
        for compr_pair in compr_pairs[:1]:
            for norm_col in norm_cols[:]:
                for genomic_region in genomic_regions[:]:
                    diff_score1,genomic_region1 = return_differential_values(indir,antibody_pair[0],compr_pair,norm_col,genomic_region)
                    diff_score2,genomic_region2 = return_differential_values(indir,antibody_pair[1],compr_pair,norm_col,genomic_region)
                    # x = diff_score1.loc[peak_overlapped_UDHS.index]
                    # y = diff_score2.loc[peak_overlapped_UDHS.index]
                    shared_index = diff_score1.index.intersection(diff_score2.index)
                    print(antibody_pair,compr_pair,len(shared_index))
                    
                    x = diff_score1.loc[shared_index][compr_pair].values
                    y = diff_score2.loc[shared_index][compr_pair].values
            
                    plt.figure(figsize=(3,3))
                    plt.scatter(x,y,rasterized=True,s=7,c='grey')
                    
                    # highlight those UTX overlapped regions                
                    x = diff_score1.loc[overlapped_ids][compr_pair].values
                    y = diff_score2.loc[overlapped_ids][compr_pair].values
                    # A
                    # data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
                    # z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
                    ## B, super slow but more accurate for cases with outliers
                    # xy = np.vstack([x,y])
                    # z = gaussian_kde(xy)(xy)
                    # z[np.where(np.isnan(z))] = 0.0
                    # idx = z.argsort()
                    # x, y, z = x[idx], y[idx], z[idx]
                    # g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
                    plt.scatter(x,y,rasterized=True,s=7,c='salmon')
                    
                    plt.axhline(y=0,c='k',lw=1,ls='--')
                    plt.axvline(x=0,c='k',lw=1,ls='--')
                    plt.xlabel('$\Delta${}  +/-{}'.format(antibody_pair[0],genomic_region1[3:]))
                    plt.ylabel('$\Delta${}  +/-{}'.format(antibody_pair[1],genomic_region2[3:]))
                    plt.title('{}\n{}\n#{}'.format(compr_pair,overlap_type,len(overlapped_ids)))
                    plt.savefig(outdir+os.sep+'{}_vs_{}_compr_{}_{}_{}_at_{}.pdf'.format(antibody_pair[0],antibody_pair[1],compr_pair,norm_col,genomic_region,overlap_type),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
                    plt.show()
                    plt.close()
                    
        # save  the info
        # shared_df = utx_peak_df.loc[shared_index]
        # shared_df.insert(3,'id',shared_df.index)
        # shared_df = shared_df.iloc[:,:4]
        # shared_df.insert(4,'x',x)
        # shared_df.insert(5,'y',y)
        # shared_df.to_csv(outdir+os.sep+'{}_{}_vs_{}.bed'.format(genomic_region,antibody_pair[0],antibody_pair[1]),sep='\t',header=None,index=False)






