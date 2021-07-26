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



def return_chipseq_rpkm_csv(celltype,antibody,genomic_region):
    # list of all the dirs with rpkm csv files
    # antibody in ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
    # genomic_region in ['Promoter','UDHS']
    file_dir1='{}/f7_chipseq/data_RPKM/1st_submission_K4me3_MLL4SC/rpkm_csv'.format(project_dir)
    file_dir2='{}/f7_chipseq/data_RPKM/202011_UTX/rpkm_csv'.format(project_dir)
    file_dir3='{}/f7_chipseq/data_RPKM/202012_K4me2/rpkm_csv'.format(project_dir)
    file_dir4='{}/f7_chipseq/data_RPKM/202102_UTX_K4me1_K27ac_K27me3/rpkm_csv'.format(project_dir)
    
    for file_dir in [file_dir1,file_dir2,file_dir3,file_dir4]:
        csv_file='{}/{}_{}_on_{}.csv'.format(file_dir,celltype,antibody,genomic_region)
        if os.path.isfile(csv_file):
            df = pd.read_csv(csv_file,sep='\t',index_col=0)
            df_rpkm = df[['RPKM']].rename(columns={"RPKM":celltype})
            df_rpkm = df_rpkm.clip(upper=np.percentile(df_rpkm,99.999))
            return df_rpkm


def return_differential_score(df,treatment,control,antibody,genomic_region,outdir,peak_overlapped_df):
    norm_t = np.sqrt(df[treatment]/df[treatment].mean())
    norm_c = np.sqrt(df[control]/df[control].mean())
    add_col='{}_over_{}'.format(treatment,control)
    df[add_col] = norm_t - norm_c
    # df = df.loc[peak_overlapped_df.index]
    
    # MA plot
    x = np.log10(df[[treatment,control]].mean(axis=1)+.1)
    y = df[add_col]
    # data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
    # z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
    # z[np.where(np.isnan(z))] = 0.0
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    plt.figure(figsize=(3,3))
    plt.scatter(x,y,rasterized=True,s=7,c='k')
    # g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=7,rasterized=True)
    plt.axhline(y=0,c='gray')
    plt.xlabel('log$_{{10}}$ (average RPKM)')
    plt.ylabel('Differential {} binding\n {} over {}'.format(antibody,treatment,control))
    plt.title(genomic_region)
    plt.savefig(outdir+os.sep+'figs_{}_{}_{}_MA.pdf'.format(antibody,genomic_region,add_col),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    return df





project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
cellTypes = ['Vector','WT','DEL','EIF']
antibodys = ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
antibodys = ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac']
genomic_regions = ['UDHS','Promoter']

outdir='f1b_differential_RPKM_MA'
os.makedirs(outdir,exist_ok=True)

# for celltype in cellTypes:
for genomic_region in genomic_regions[:]:
    for antibody in antibodys[:]:
        # only keep those region overlap with either peak of chip-seq data
        peak_overlapped_file = '{}/f7_chipseq/data_peak_overlap/overlapped/{}_overlap_{}.bed'.format(project_dir,genomic_region,antibody)
        peak_overlapped_df = pd.read_csv(peak_overlapped_file,sep='\t',index_col=3,header=None)
        df_vector = return_chipseq_rpkm_csv('Vector',antibody,genomic_region)
        df_wt = return_chipseq_rpkm_csv('WT',antibody,genomic_region)
        df_del = return_chipseq_rpkm_csv('DEL',antibody,genomic_region)
        df_eif = return_chipseq_rpkm_csv('EIF',antibody,genomic_region)
        df = pd.concat([df_vector,df_wt,df_del,df_eif],axis=1).loc[peak_overlapped_df.index]
        print(genomic_region,antibody,df.shape)
        df = return_differential_score(df,'WT','Vector',antibody,genomic_region,outdir,peak_overlapped_df)
        df = return_differential_score(df,'DEL','WT',antibody,genomic_region,outdir,peak_overlapped_df)
        df = return_differential_score(df,'EIF','DEL',antibody,genomic_region,outdir,peak_overlapped_df)
#         print(celltype,antibody,df);exit()
        df.to_csv(outdir+os.sep+'{}_{}_differential_RPKM.csv'.format(antibody,genomic_region))







