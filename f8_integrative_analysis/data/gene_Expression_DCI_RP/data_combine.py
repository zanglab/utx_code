import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats

    
## ==== main
project_dir = "/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir = "/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
factors = ['UTX','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']


# ==== combine the data
outdir = 'fz_data_combined'
os.makedirs(outdir,exist_ok=True)

# read the expression/DCI data
expression_DCI_file = '{}/f8_integrative_analysis/data/gene_Expression_DCI/fz_data_combined/summary_TPM_DESeq2_DCI.csv'.format(project_dir)
expression_DCI_df = pd.read_csv(expression_DCI_file,index_col=0)

# dropna_kept_col = [i for i in expression_DCI_df.columns if not re.search('padj',i)]
# dropna_kept_index = expression_DCI_df[dropna_kept_col].dropna().index
# expression_DCI_df = expression_DCI_df.loc[dropna_kept_index]

# read the HM RPKM/RP data
rpkm_rp_file = '{}/f7_chipseq/f9_diff_binding_on_promoters/f0_integrate_RPKM_RP/f1_combined_data/combined_RPKM_RP_on_Promoter.csv'.format(project_dir)
rpkm_rp_df = pd.read_csv(rpkm_rp_file,index_col=0)

# combine the data
combined_df = pd.concat([expression_DCI_df,rpkm_rp_df],axis = 1)
combined_df.to_csv(outdir+os.sep+'TPM_DEseq2_DCI_RPKM_RP.csv')


# write out the RP of k4me1, k4me3, and k27ac
rp_cols = [i for i in rpkm_rp_df.columns if not i.endswith('log2AVG')]
rp_cols = [i for i in rp_cols if i.startswith('RP_')]
rp_cols = [i for i in rp_cols if re.search('H3K4me1|H3K4me3|H3K27ac',i)]
rp_cols = np.append(expression_DCI_df.columns,rp_cols)
writer = pd.ExcelWriter(outdir+os.sep+'TPM_DEseq2_DCI_RPKM_RP.xlsx')
combined_df[rp_cols].to_excel(writer,'All genes')
writer.save()


writer = pd.ExcelWriter(outdir+os.sep+'TPM_DEseq2_DCI.xlsx')
expression_DCI_df.to_excel(writer,'All genes')
writer.save()




