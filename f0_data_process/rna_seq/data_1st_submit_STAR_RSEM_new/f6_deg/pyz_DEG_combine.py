import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats


deseq2_labels = {'WT_over_Vector':'treated_WT_vs_ctrl_Vector.deseq2.csv',\
                 'DEL_over_WT':'treated_del_cIDR_vs_ctrl_WT.deseq2.csv',\
                 'EIF_over_DEL':'treated_UTX_eIFIDR_vs_ctrl_del_cIDR.deseq2.csv',\
                 'TPR_over_WT':'treated_del_TPR_vs_ctrl_WT.deseq2.csv'}

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)

outdir = 'fz_deseq2_out_combined'
os.makedirs(outdir,exist_ok=True)

compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL','TPR_over_WT']
combined_df = pd.DataFrame()

for compr_type in compr_types[:]:
    deg_file = expr_dir+os.sep+deseq2_labels[compr_type]
    deg_df = pd.read_csv(deg_file,index_col=0)
    deg_tmp = deg_df[['log2FoldChange','padj']]
    deg_tmp.log2FoldChange = deg_tmp.log2FoldChange.round(6)
    deg_tmp.columns = ['{}_log2FoldChange'.format(compr_type),'{}_padj'.format(compr_type)]
    combined_df = pd.concat([combined_df,deg_tmp],axis=1)
combined_df.to_csv(outdir+os.sep+'deseq2_combined.csv')
