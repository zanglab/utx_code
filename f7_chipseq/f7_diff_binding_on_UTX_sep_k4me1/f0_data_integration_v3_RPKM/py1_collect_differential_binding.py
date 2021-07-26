import sys,argparse
import os,glob,re
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



def return_chipseq_rpkm_csv(celltype,peak_file,factor,genomic_region,norm_col):
    # list of all the dirs with rpkm csv files
    # factor in ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
    file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_RPKM/rpkm_csv'.format(project_dir)
    csv_file='{}/{}_{}_on_{}_{}.csv'.format(file_dir,celltype,factor,peak_file,genomic_region)
    print(celltype,peak_file,factor,genomic_region,norm_col)
    df = pd.read_csv(csv_file,sep='\t',index_col=0) 
    df_rpkm = df[['ReadCount']].rename(columns={"ReadCount":'{}_{}'.format(factor,celltype)})
    df_rpkm = df_rpkm.clip(upper=np.percentile(df_rpkm,99.99))
    
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    # norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),'total']/1000000
    # norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),'total_in_islads']/1000000
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    if genomic_region=='es2kb':
        norm_factor = norm_factor*4 # normalized reads per kb
    return df_rpkm/norm_factor


def return_differential_score(df,treatment,control,factor,genomic_region,outdir):
    # == log2FC
    pc=0.1
    # add_col='{}_over_{}_log2FC'.format(treatment,control)
    # df[add_col] = np.log2((df[treatment]+pc)/(df[control]+pc))
    # MA plot
    # x = df[[treatment,control]].mean(axis=1)
    treatment = '{}_{}'.format(factor,treatment)
    control = '{}_{}'.format(factor,control)
    x = np.log2(df[[treatment,control]].mean(axis=1)+pc)
    y = np.log2((df[treatment]+pc)/(df[control]+pc))
    df['{}_over_{}_log2AVG'.format(treatment,control)] = x
    df['{}_over_{}_log2FC'.format(treatment,control)] = y
    
    return df





project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# cellTypes = ['Vector','WT','DEL','EIF','MT2','TPR']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
# genomic_regions = ['es500bp','es2kb',]
# peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
# norm_cols = ['total','total_in_islads']
genomic_regions = ['es2kb',]
peak_files = ['UTX_peaks']
norm_cols = ['total']

outdir='f1_differential_NormReadCount'
os.makedirs(outdir,exist_ok=True)

# for celltype in cellTypes:
for peak_file in peak_files[:]:
    for factor in factors[:]:
        for genomic_region in genomic_regions[:]:
            for norm_col in norm_cols[:]:
                df_vector = return_chipseq_rpkm_csv('Vector',peak_file,factor,genomic_region,norm_col)
                df_wt = return_chipseq_rpkm_csv('WT',peak_file,factor,genomic_region,norm_col)
                df_del = return_chipseq_rpkm_csv('DEL',peak_file,factor,genomic_region,norm_col)
                # test if MT2 file exist
                file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_RPKM/rpkm_csv'.format(project_dir)
                test_file='{}/{}_{}_on_{}_{}.csv'.format(file_dir,'MT2',factor,peak_file,genomic_region)
                if os.path.isfile(test_file):
                    df_mt2 = return_chipseq_rpkm_csv('MT2',peak_file,factor,genomic_region,norm_col)
                    df = pd.concat([df_vector,df_wt,df_del,df_mt2],axis=1)
                    df = return_differential_score(df,'WT','Vector',factor,genomic_region,outdir)
                    df = return_differential_score(df,'DEL','WT',factor,genomic_region,outdir)
                    df = return_differential_score(df,'MT2','DEL',factor,genomic_region,outdir)
                else:
                    df = pd.concat([df_vector,df_wt,df_del],axis=1)
                    df = return_differential_score(df,'WT','Vector',factor,genomic_region,outdir)
                    df = return_differential_score(df,'DEL','WT',factor,genomic_region,outdir)
                # print(celltype,factor,df);exit()
                # for SICER-df islands, keep only those UTX islands with increased signal in WT than in Vector
                if re.search('islands',peak_file):
                    islands_increased_file='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/{}_increased.bed'.format(project_dir,peak_file)
                    islands_increased_df = pd.read_csv(islands_increased_file,sep='\t',header=None)
                    kept_ids = ['_'.join(ii) for ii in islands_increased_df.loc[:,:2].astype(str).values]
                    df.loc[kept_ids].round(8).to_csv(outdir+os.sep+'{}_{}_DiffNormReads_on_{}_by_{}.csv'.format(genomic_region,factor,peak_file,norm_col))
                else:
                    df.round(8).to_csv(outdir+os.sep+'{}_{}_DiffNormReads_on_{}_by_{}.csv'.format(genomic_region,factor,peak_file,norm_col))







