import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect




# ==== main 
    
outdir = 'f1_RPKM_collection'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'

sub_dir_names=['data_1st_submit','data_20201127_1205_merged','data_20201209','data_20201220']
binding_types=['Promoter','Promoter_es10kb']

for binding_type in binding_types:
    # init the RPKM df for each binding_type
    df_rpkm = pd.DataFrame()   
    for sub_dir_name in sub_dir_names:
        sub_dir = '{}/f3_differential_chip_binding/{}/chip_binding_RPKM_promoter_UDHS/rpkm_csv'.format(project_dir,sub_dir_name)
        csv_files = glob.glob('{}/*{}.csv'.format(sub_dir,binding_type))
        for csv_file in csv_files:
            basename = os.path.basename(csv_file).split('_on_{}'.format(binding_type))[0]
            # get the RPKM and read count
            df = pd.read_csv(csv_file,sep='\t',index_col=0)
            df_rpkm['{}_{}'.format(sub_dir_name,basename)] = df['RPKM']
            # print(df_rpkm);exit(0)
            print('{}_{}'.format(sub_dir_name,basename))

    df_rpkm.to_csv(outdir+os.sep+'{}_RPKM.csv'.format(binding_type))



    
    
