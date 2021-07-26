import sys,argparse
import os,glob
import numpy as np
import pandas as pd




project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

peak_overlapped_UDHS='{}/f7_chipseq/data_peak_overlap/overlapped/UDHS_overlap_ALL.bed'.format(project_dir)
peak_overlapped_UDHS_df = pd.read_csv(peak_overlapped_UDHS,sep='\t',header=None,index_col=3)
# print(peak_overlapped_UDHS_df);exit()

outdir='rpkm_csv'
os.makedirs(outdir,exist_ok=True)

file_dirs=['1st_submission_K4me3_MLL4SC','202011_UTX','202012_K4me2','202102_UTX_K4me1_K27ac_K27me3']
file_dirs=['1st_submission_K4me3_MLL4SC','202011_UTX','202012_K4me2']

for file_dir in file_dirs:
    outdir_tmp = '{}/{}'.format(outdir,file_dir)
    os.makedirs(outdir_tmp,exist_ok=True)
    
    csv_files = glob.glob('../{}/rpkm_csv/*UDHS.csv'.format(file_dir))
    for csv_file in csv_files:
        outname = os.path.basename(csv_file)
        df = pd.read_csv(csv_file,sep='\t',index_col=0)
        df = df.loc[peak_overlapped_UDHS_df.index]
        df.to_csv('{}/{}'.format(outdir_tmp,outname),sep='\t')
        print(outdir_tmp,outname)








