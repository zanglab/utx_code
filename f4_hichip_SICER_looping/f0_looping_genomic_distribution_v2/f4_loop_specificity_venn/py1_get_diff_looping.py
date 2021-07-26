import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

    

def get_diff_looping(hm,cellTypes, gene_info_columns,indir,outdir,subdir):

    os.makedirs(outdir+os.sep+subdir,exist_ok=True)
    df = pd.DataFrame() 
    for celltype in cellTypes:
        hm_df = pd.read_csv('{}/{}/{}_{}.csv'.format(indir,subdir,hm,celltype),index_col=0)
        hm_df['cell_type']=celltype
        hm_df['cell_type_count']=1
        df = pd.concat([df,hm_df])
    
    df_cell = df.groupby(df.index)['cell_type'].apply(','.join)    
    df_cell_count = df.groupby(df.index)['cell_type_count'].apply(sum)    
    df_overlap_info = df.groupby(df.index)[gene_info_columns].first()
    pd.concat([df_cell,df_cell_count,df_overlap_info],axis=1).to_csv(outdir+os.sep+subdir+os.sep+'{}.csv'.format(hm))      



# ==== main()


indir = '../f2_loop_anchor_annotation/f2_loop_genomic_feature'
outdir = 'f1_loop_specificity'
subdirs=['data_1st_submission_rep_combined','data_202008']

gene_info_columns = ['anchor1_promoter_Overlapped_IDs', 
                      'anchor2_promoter_Overlapped_IDs']

# plot the figs
subdir=subdirs[0]
hm='H3K4me3'
cellTypes = ['Vector','WT','DEL','EIF']   
get_diff_looping(hm,cellTypes, gene_info_columns,indir,outdir,subdir)
   

# ==== batch2 H3K4me3
subdir=subdirs[1]
hm='H3K4me3'
cellTypes = ['Vector','WT','DEL','EIF','TPR']    
get_diff_looping(hm,cellTypes, gene_info_columns,indir,outdir,subdir)

# ==== batch2 H3K27ac
hm='H3K27ac'
get_diff_looping(hm,cellTypes, gene_info_columns,indir,outdir,subdir)


