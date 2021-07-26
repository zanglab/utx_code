import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

    
indir = 'f2_loop_genomic_feature'
outdir = 'f3_loop_genomic_distribution_data'
subdirs=['data_1st_submission_rep_combined','data_1st_submission_sep_rep','data_202008']

for subdir in subdirs[:]:
    os.makedirs('{}/{}'.format(outdir,subdir),exist_ok=True)
    all_files = glob.glob('{}/{}/*'.format(indir,subdir))
    prenames = np.unique([os.path.basename(i).split('.csv')[0] for i in all_files])
    
    # colloect the genomic features of anchors
    sum_df = pd.DataFrame()
    for hm in prenames[:]:
        df = pd.read_csv('{}/{}/{}.csv'.format(indir,subdir,hm),index_col=0)
        
        # promoter vs. enenhancer (simplified as none gene body region)
        PP = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
        PE = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_UDHS_IfOverlap']==1)].shape[0]
        PO = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_UDHS_IfOverlap']!=1)].shape[0]
        EP = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_UDHS_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
        OP = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_UDHS_IfOverlap']!=1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
        EE = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_UDHS_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_UDHS_IfOverlap']==1)].shape[0]
        
        # summarize the interact information
        sum_df.loc['PP',hm] = PP
        sum_df.loc['PE',hm] = PE
        sum_df.loc['PO',hm] = PO
        sum_df.loc['EP',hm] = EP
        sum_df.loc['OP',hm] = OP
        sum_df.loc['EE',hm] = EE
        
        sum_df.loc['P-P',hm] = PP
        sum_df.loc['P-E',hm] = PE+EP
        sum_df.loc['P-O',hm] = PO+OP
        sum_df.loc['E-E',hm] = EE
        
        sum_df.loc['Total',hm] = df.shape[0]
        sum_df['{}%'.format(hm)] = sum_df[hm]/sum_df.at['Total',hm]
    
    sum_df.to_csv(outdir+os.sep+subdir+os.sep+'loop_distribution_summary.csv')    
    

