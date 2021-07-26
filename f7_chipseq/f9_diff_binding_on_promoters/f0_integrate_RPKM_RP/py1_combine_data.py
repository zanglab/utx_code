import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=18
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks")
# from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde



def return_chipseq_rpkm_csv(celltype,factor,count_type,genomic_region,col_name):
    # list of all the dirs with rpkm csv files
    # factor in ['UTX','MLL4','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3']
    file_dir='{}/f7_chipseq/f9_diff_binding_on_promoters/data_{}/{}_csv'.format(project_dir,count_type,count_type)
    csv_file='{}/{}_{}_on_{}.csv'.format(file_dir,celltype,factor,genomic_region)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = df[[col_name]].rename(columns={col_name:'{}_{}_{}'.format(count_type,factor,celltype)})
    df = df.clip(upper=np.percentile(df,99.99))
    return df


def return_differential_score(df,treatment,control,factor,count_type):
    # == log2FC
    pc=0.1
    t = '{}_{}_{}'.format(count_type,factor,treatment)
    c = '{}_{}_{}'.format(count_type,factor,control)
    x = np.log2(df[[t,c]].mean(axis=1)+pc)
    y = np.log2((df[t]+pc)/(df[c]+pc))
    df['{}_{}_{}_over_{}_log2AVG'.format(count_type,factor,treatment,control)] = x
    df['{}_{}_{}_over_{}_log2FC'.format(count_type,factor,treatment,control)] = y
    return df




# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
count_type_matches = {'RPKM':['Promoter_es2kb','RPKM'],
                      'RP':['Promoter','RP_M']}

outdir='f1_combined_data'
os.makedirs(outdir,exist_ok=True)


df_combined = pd.DataFrame()
for factor in factors[:]:
    ## get the RPKM for each factor
    for count_type in count_type_matches.keys():
        genomic_region = count_type_matches[count_type][0]
        col_name = count_type_matches[count_type][1]
        
        df_vector = return_chipseq_rpkm_csv('Vector',factor,count_type,genomic_region,col_name)
        df_wt = return_chipseq_rpkm_csv('WT',factor,count_type,genomic_region,col_name)
        df_del = return_chipseq_rpkm_csv('DEL',factor,count_type,genomic_region,col_name)
        df_eif = return_chipseq_rpkm_csv('EIF',factor,count_type,genomic_region,col_name)
        
        df_factor = pd.concat([df_vector,df_wt,df_del,df_eif],axis=1)
        df_factor = return_differential_score(df_factor,'WT','Vector',factor,count_type)
        df_factor = return_differential_score(df_factor,'DEL','WT',factor,count_type)
        df_factor = return_differential_score(df_factor,'EIF','DEL',factor,count_type)
        df_combined = pd.concat([df_combined,df_factor],axis=1)

df_combined.round(6).to_csv(outdir+os.sep+'combined_RPKM_RP_on_Promoter.csv')



# x='RPKM_H3K27ac_WT_over_Vector_log2FC'
# y='RP_H3K27ac_WT_over_Vector_log2FC'
# plt.scatter(df_combined[x],df_combined[y])



