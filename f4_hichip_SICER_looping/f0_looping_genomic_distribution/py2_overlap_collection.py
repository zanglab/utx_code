import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"

    

outdir = 'f2_loop_genomic_feature'
os.makedirs(outdir,exist_ok=True)

hms=['K27ACDEL','K27ACEIF','K27ACTPR','K27ACVEC','K27ACWT','K4M3DEL3','K4M3DTPR','K4M3EIF','K4M3PCDH','K4M3WT']
for hm in hms:
    hm_df = pd.read_csv('data_loop_bedpe/{}.bedpe'.format(hm),sep='\t')
    hm_df.index = hm_df.index+1 # index starts from 1
    for anchor_flag in ['anchor1','anchor2']:
        for genomic_feature in ['promoter','exon','intron','genebody']:
            # read the overlapped files
            anchor_file='overlapped/{}_{}_{}.bed'.format(hm,genomic_feature,anchor_flag)
            anchor_df = pd.read_csv(anchor_file,sep='\t',index_col=4)
            anchor_df = anchor_df[['IfOverlap','Overlapped_genes']]
            col_names = ['{}_{}_IfOverlap'.format(anchor_flag,genomic_feature),'{}_{}_Overlapped_genes'.format(anchor_flag,genomic_feature)]    
            anchor_df.columns = col_names
            hm_df = pd.concat([hm_df,anchor_df],axis=1)
    hm_df.to_csv(outdir+os.sep+'{}.csv'.format(hm))

