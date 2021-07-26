import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

    

def get_diff_looping(hms,cell_types, gene_info_columns,indir, outdir,flag):

    columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    df = pd.DataFrame()
    for ii in np.arange(len(hms)):
        hm = hms[ii]
        cell_type = cell_types[ii]
        hm_df = pd.read_csv('{}/{}.csv'.format(indir,hm),index_col=0)
        # keep only those high-confident loops
        hm_df = hm_df[hm_df['fdr']<1e-5]
        hm_df['cell_type']=cell_type
        hm_df['cell_type_count']=1
        df = pd.concat([df,hm_df])

    df_cell = df.groupby(columns)['cell_type'].apply(','.join)    
    df_cell_count = df.groupby(columns)['cell_type_count'].apply(sum)    
    df_overlap_info = df.groupby(columns)[gene_info_columns].first()
    pd.concat([df_cell,df_cell_count,df_overlap_info],axis=1).to_csv(outdir+os.sep+'{}.csv'.format(flag))      




# ==== main()

indir = '../f0_looping_genomic_distribution/f2_loop_genomic_feature_new_with_UDHS'
outdir = 'f2_diff_looping_high_thre'
os.makedirs(outdir,exist_ok=True)

gene_info_columns = ['anchor1_promoter_Overlapped_genes', 'anchor1_udhs_IfOverlap', 
                     'anchor2_promoter_Overlapped_genes', 'anchor2_udhs_IfOverlap',]

cell_types =['VEC','WT','DEL','EIF','TPR']

hms = ['K27ACVEC','K27ACWT','K27ACDEL','K27ACEIF','K27ACTPR',]    
get_diff_looping(hms,cell_types,gene_info_columns, indir, outdir, 'H3K27ac')
    
hms = ['K4M3PCDH','K4M3WT','K4M3DEL3','K4M3EIF','K4M3DTPR',]
get_diff_looping(hms,cell_types,gene_info_columns, indir, outdir, 'H3K4me3')

