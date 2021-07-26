import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
    
    

# ==== main() 

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

    
# project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)

outdir = 'fz_utx_binding_DCI_combined'
os.makedirs(outdir,exist_ok=True)


peak_files = ['UTX_peaks','UTX_islands','UTXFEB_islands','UTXFEB_peaks']
subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']

name_match_df = pd.read_excel('dci_utx_matching.xlsx',index_col=0)


for peak_file in peak_files[:]:
    # for each hichip data type and each utx binding file
    combined_df = pd.DataFrame()
    for subdir in subdirs[:]:
        for dci_file_basename in name_match_df.index[:]:
            dci_file = 'f1_UTX_binding_DCI/{}/{}_{}_DCI.csv'.format(subdir,dci_file_basename,peak_file)
            if os.path.isfile(dci_file):
                # print(subdir,dci_file_basename,peak_file)
                utx_df = pd.read_csv(dci_file,sep='\t',index_col=4)
                compr_figname = name_match_df.loc[dci_file_basename,'compr_figname']
                hm = name_match_df.loc[dci_file_basename,'hm']
                col_rename = '{}--{}--{}'.format(subdir,hm,compr_figname)
                df = utx_df[['info']].rename(columns={'info':col_rename})
                combined_df = pd.concat([combined_df,df],axis=1)
    combined_df.to_csv('{}/{}_HiChIP_DCI.csv'.format(outdir,peak_file))
    # hichip_data, hm, compr_celltypes = column.split('--')






