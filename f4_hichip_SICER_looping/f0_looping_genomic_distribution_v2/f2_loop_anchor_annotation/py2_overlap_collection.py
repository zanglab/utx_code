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


indir = 'f1_overlapped'
outdir = 'f2_loop_genomic_feature'
subdirs=['data_1st_submission_rep_combined','data_1st_submission_sep_rep','data_202008']

for subdir in subdirs[:]:
    os.makedirs('{}/{}'.format(outdir,subdir),exist_ok=True)
    all_files = glob.glob('{}/{}/*'.format(indir,subdir))
    prenames = np.unique([os.path.basename(i).split('_anchor')[0] for i in all_files])
    # for each type of data, merge the anchor1 and anchor2 files
    for prename in prenames[:]:
        hm_df = pd.DataFrame()
        for genomic_feature in ['promoter','UDHS'][:]:
            for anchor in ['anchor1','anchor2']:
                df_file = '{}/{}//{}_{}_{}.bed'.format(indir,subdir,prename,anchor,genomic_feature)
                df = pd.read_csv(df_file,sep='\t',index_col=4)
                df = df[['IfOverlap', 'Overlapped_genes']]
                col_names = ['{}_{}_IfOverlap'.format(anchor,genomic_feature),'{}_{}_Overlapped_IDs'.format(anchor,genomic_feature)]    
                df.columns = col_names
                hm_df = pd.concat([hm_df,df],axis=1)
            
        hm_df.to_csv('{}/{}/{}.csv'.format(outdir,subdir,prename))




