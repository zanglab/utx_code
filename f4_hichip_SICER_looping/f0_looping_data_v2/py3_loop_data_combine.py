import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats



def add_coordinates(df):
    coordinates = pd.DataFrame([i.split('_') for i in df.index])
    coordinates.index = df.index
    coordinates.columns = ['chr1','start1','chr2','start2']
    coordinates['start1'] = coordinates.start1.astype(int)
    coordinates['end1'] = coordinates.start1+5000
    coordinates['start2'] = coordinates.start2.astype(int)
    coordinates['end2'] = coordinates.start2+5000
    df = pd.concat([coordinates[['chr1','start1','end1','chr2','start2','end2']],df],axis=1)
    return df



writer = pd.ExcelWriter('summary_HiChIP_loops.xlsx')
norm_col = 'number_of_pairs_after_duplicate_removal'

# ==== batch1 H3K4me3
norm_file = '/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_1st_submit_sicer2_merged_islands_new/qc_summary_rename.xlsx'
norm_df = pd.read_excel(norm_file,index_col=1)
indir='data_1st_submission_reindex'
outname='batch1_H3K4me3'
hm='H3K4me3'
cell_types =['Vector','WT','DEL','EIF']
df_combined=pd.DataFrame()
for celltype in cell_types[:]:
    for rep in ['rep1','rep2']:
        df = pd.read_csv('{}/{}_{}_{}.bedpe'.format(indir,hm,celltype,rep),index_col=0,sep='\t')
        norm_factor = norm_df.loc[norm_col,'{}_{}_{}'.format(hm,celltype,rep)]/100000000
        df['normalized_count'] = (df['count'].values/norm_factor).round(2)
        df = df[['count','normalized_count','fdr']]
        df.columns = ['{} {} raw count'.format(celltype,rep),'{} {} normalized count'.format(celltype,rep),'{} {} fdr'.format(celltype,rep)]
        # df = dfs[['normalized_count']].rename(columns={'normalized_count':'{} {}'.format(celltype,rep)})
        df_combined = pd.concat([df_combined,df],axis=1)
# save to excel
df_combined = add_coordinates(df_combined)
df_combined.to_excel(writer,outname,index=False)



# ==== batch2 H3K4me3
norm_file = '/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_202008_sicer2_merged_islands_with_new_k27/qc_summary_rename.xlsx'
norm_df = pd.read_excel(norm_file,index_col=1)
indir='data_202008_reindex'
outname='batch2_H3K4me3'
hm='H3K4me3'
cell_types =['Vector','WT','DEL','EIF','TPR']
df_combined=pd.DataFrame()
for celltype in cell_types[:]:
    df = pd.read_csv('{}/{}_{}.bedpe'.format(indir,hm,celltype),index_col=0,sep='\t')
    norm_factor = norm_df.loc[norm_col,'{}_{}'.format(hm,celltype)]/100000000
    df['normalized_count'] = (df['count'].values/norm_factor).round(2)
    df = df[['count','normalized_count','fdr']]
    df.columns = ['{} raw count'.format(celltype),'{} normalized count'.format(celltype),'{} fdr'.format(celltype)]
    # df = df[['normalized_count']].rename(columns={'normalized_count':'{}'.format(celltype)})
    df_combined = pd.concat([df_combined,df],axis=1)
# save to excel
df_combined = add_coordinates(df_combined)
df_combined.to_excel(writer,outname,index=False)



# ==== batch2 H3K27ac
indir='data_202008_reindex'
outname='batch2_H3K27ac'
hm='H3K27ac'
cell_types =['Vector','WT','DEL','EIF','TPR']
df_combined=pd.DataFrame()
for celltype in cell_types[:]:
    df = pd.read_csv('{}/{}_{}.bedpe'.format(indir,hm,celltype),index_col=0,sep='\t')
    norm_factor = norm_df.loc[norm_col,'{}_{}'.format(hm,celltype)]/100000000
    df['normalized_count'] = (df['count'].values/norm_factor).round(2)
    df = df[['count','normalized_count','fdr']]
    df.columns = ['{} raw count'.format(celltype),'{} normalized count'.format(celltype),'{} fdr'.format(celltype)]
    # df = df[['normalized_count']].rename(columns={'normalized_count':'{}'.format(celltype)})
    df_combined = pd.concat([df_combined,df],axis=1)
# save to excel
df_combined = add_coordinates(df_combined)
df_combined.to_excel(writer,outname,index=False)


writer.save()


