import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect


#for indir in ['data_1st_submission','data_202008']:
for indir in ['data_202008_new']:
    outdir='{}_reindex'.format(indir)
    os.makedirs(outdir,exist_ok=True)
    bedpe_files = glob.glob('{}/*.bedpe'.format(indir))
    for bedpe_file in bedpe_files:
        basename = os.path.basename(bedpe_file)
        df = pd.read_csv(bedpe_file,sep='\t')
        pos1_index = df['chr1']+'_'+df['start1'].astype(str)
        pos2_index = df['chr2']+'_'+df['start2'].astype(str)
        re_index = pos1_index+'_'+pos2_index
        df.insert(0,'pos2_index',pos2_index)
        df.insert(0,'pos1_index',pos1_index)
        df.insert(0,'re_index',re_index)
        df.to_csv('{}/{}'.format(outdir,basename),index=False,sep='\t')
        print(basename)