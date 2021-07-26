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





def differential_binding_MA(df,treatment,control,factor,peak_file,outdir):


    treatment_col = '{}_{}'.format(factor,treatment)
    control_col = '{}_{}'.format(factor,control)
    x = df['{}_over_{}_log2AVG'.format(treatment_col,control_col)] 
    y = df['{}_over_{}_log2FC'.format(treatment_col,control_col)]
    data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
    z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
    z[np.where(np.isnan(z))] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    plt.figure(figsize=(3,3))
#     plt.scatter(x,y,rasterized=True,s=7,c='k')
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=7,rasterized=True)
    plt.axhline(y=0,c='gray')
    plt.xlabel('log$_{{2}}$ (average {} binding)'.format(factor))
    plt.ylabel('Differential {} binding \n {} over {}'.format(factor,treatment,control))
    # plt.title('{} {}'.format(factor,genomic_region))
    plt.savefig(outdir+os.sep+'{}_{}_{}_over_{}_MA.pdf'.format(peak_file,factor,treatment,control,),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    





# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# cellTypes = ['Vector','WT','DEL','EIF']
# genomic_regions = ['es500bp','es2kb',]
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']

indir='../f0_data_integration/f2_combined_data/'
outdir='f1_differential_NormReads_MA'
os.makedirs(outdir,exist_ok=True)

# for celltype in cellTypes:
for peak_file in peak_files[:]:
    df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
    for factor in factors[:]:
        differential_binding_MA(df,'WT','Vector',factor,peak_file,outdir)
        # differential_binding_MA(df,'DEL','WT',factor,genomic_region,outdir)
        # differential_binding_MA(df,'EIF','DEL',factor,genomic_region,outdir)







