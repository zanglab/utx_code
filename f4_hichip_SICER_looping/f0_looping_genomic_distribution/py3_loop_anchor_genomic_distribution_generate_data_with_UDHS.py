import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

    
indir = 'f2_loop_genomic_feature_new_with_UDHS'
outdir = 'f3_loop_genomic_distribution_new_with_UDHS'
os.makedirs(outdir,exist_ok=True)

# columns = ['anchor1_promoter_IfOverlap',
#  'anchor1_exon_IfOverlap',
#  'anchor1_intron_IfOverlap',
#  'anchor1_genebody_IfOverlap',
#  'anchor2_promoter_IfOverlap',
#  'anchor2_exon_IfOverlap',
#  'anchor2_intron_IfOverlap',
#  'anchor2_genebody_IfOverlap']

# P-promoter
# E,IG - enhancer, intergenic
# GB - gene body


hms=['K27ACDEL','K27ACEIF','K27ACTPR','K27ACVEC','K27ACWT','K4M3DEL3','K4M3DTPR','K4M3EIF','K4M3PCDH','K4M3WT']
# hms=['K27ACDEL','K27ACEIF',]
sum_df = pd.DataFrame()

for hm in hms:
    df = pd.read_csv('{}/{}.csv'.format(indir,hm),index_col=0)
    # promoter vs. enenhancer (simplified as none gene body region)
    PP = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
    PE = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_udhs_IfOverlap']==1)].shape[0]
    PO = df[(df['anchor1_promoter_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_udhs_IfOverlap']!=1)].shape[0]
    EP = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_udhs_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
    OP = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_udhs_IfOverlap']!=1)&(df['anchor2_promoter_IfOverlap']==1)].shape[0]
    EE = df[(df['anchor1_promoter_IfOverlap']!=1)&(df['anchor1_udhs_IfOverlap']==1)&(df['anchor2_promoter_IfOverlap']!=1)&(df['anchor2_udhs_IfOverlap']==1)].shape[0]
    
    
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

print(sum_df)
sum_df.to_csv(outdir+os.sep+'loop_distribution_summary.csv')


