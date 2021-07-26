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



def return_chipseq_rpkm_csv(celltype,peak_file,factor,genomic_region):
    # list of all the dirs with rpkm csv files
    # factor in ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
    file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_RPKM/rpkm_csv'.format(project_dir)
    csv_file='{}/{}_{}_on_{}_{}.csv'.format(file_dir,celltype,factor,peak_file,genomic_region)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df_rpkm = df[['ReadCount']].rename(columns={"ReadCount":'{}_{}'.format(factor,celltype)})
    df_rpkm = df_rpkm.clip(upper=np.percentile(df_rpkm,99.99))
    
    return df_rpkm['{}_{}'.format(factor,celltype)]




    
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
cellTypes = ['Vector','WT','DEL','EIF']
peak_files = ['UTX_all_islands_merge','UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
genomic_region = 'es500bp'

outdir='f1_UTX_UTXFEB_binding_scatter'
os.makedirs(outdir,exist_ok=True)



for peak_file in peak_files:
    # for each pair of cell types, compr diff binding between two factors
    for celltype in cellTypes:

        utx_signal = return_chipseq_rpkm_csv(celltype,peak_file,'UTX',genomic_region)
        utxfeb_signal = return_chipseq_rpkm_csv(celltype,peak_file,'UTXFEB',genomic_region)

        x = np.log2(utx_signal+1)
        y = np.log2(utxfeb_signal+1)

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
        plt.xlabel('log2 (UTX ReadCount)')
        plt.ylabel('log2 (UTXFEB ReadCount)')
        plt.title('{}'.format(celltype))
        plt.savefig(outdir+os.sep+'{}_UTX_vs_UTXFEB_in_{}.pdf'.format(peak_file,celltype),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
            
    





