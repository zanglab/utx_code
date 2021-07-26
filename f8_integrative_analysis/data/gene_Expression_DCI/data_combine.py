import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats


outdir = 'fz_data_combined'
os.makedirs(outdir,exist_ok=True)


project_dir = "/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir = "/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
deseq2_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
tpm_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f3_expr'.format(project_dir)
dci_1st_dir='{}/f5_hichip/f1_hichip_bart3d_new/f1_DEG_promoter_DCI/f1_promoter_DCI/bart3d_dis200k_data_1st_submit'.format(project_dir)
dci_202008_dir='{}/f5_hichip/f1_hichip_bart3d_new/f1_DEG_promoter_DCI/f1_promoter_DCI/bart3d_dis200k_data202008'.format(project_dir)


# ==== read the Deseq2 results
deseq2_df = pd.read_csv('{}/deseq2_combined.csv'.format(deseq2_dir),index_col=0)


# ==== read the gene expression TPM data
tpm_df = pd.read_csv('{}/tpm.csv'.format(tpm_dir),index_col=0)
col_rename = {'UTX_FUSIDR':'FUS',
          'UTX_eIFIDR':'EIF',
          'del_TPR':'TPR',
          'del_cIDR':'DEL'}
new_columns = []
for col in tpm_df.columns:
    celltype,rep = col.split('_rep')
    if celltype in col_rename.keys():
        new_columns.append('{}_rep{}_TPM'.format(col_rename[celltype],rep))
    else:
        new_columns.append('{}_TPM'.format(col))
tpm_df.columns = new_columns
tpm_df = tpm_df[['Vector_rep1_TPM', 'Vector_rep2_TPM', 
                'WT_rep1_TPM', 'WT_rep2_TPM',
                'DEL_rep1_TPM', 'DEL_rep2_TPM',
                'EIF_rep1_TPM', 'EIF_rep2_TPM',
                'TPR_rep1_TPM','TPR_rep2_TPM']]


# ==== read the DCI profile from batch1 data
compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']
hms = ['H3K4me3']
combined_dci_batch1= pd.DataFrame()
for hm in hms:
    for compr_type in compr_types:
        dci_prename = '{}_{}'.format(hm,compr_type)
        dci_df = pd.read_csv('{}/{}_promoter_DCI.csv'.format(dci_1st_dir,dci_prename),index_col=4,sep='\t')
        dci_df = dci_df[['info']].rename(columns={'info':'Batch1_{}_DCI'.format(dci_prename)}).round(6)
        combined_dci_batch1 = pd.concat([combined_dci_batch1,dci_df],axis=1)
        


# ==== read the DCI profile from batch2 data
compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL','TPR_over_WT']
hms = ['H3K4me3','H3K27ac']
combined_dci_batch2= pd.DataFrame()
for hm in hms:
    for compr_type in compr_types:
        dci_prename = '{}_{}'.format(hm,compr_type)
        dci_df = pd.read_csv('{}/{}_promoter_DCI.csv'.format(dci_202008_dir,dci_prename),index_col=4,sep='\t')
        dci_df = dci_df[['info']].rename(columns={'info':'Batch2_{}_DCI'.format(dci_prename)}).round(6)
        combined_dci_batch2 = pd.concat([combined_dci_batch2,dci_df],axis=1)
        



# create the Excel file
writer = pd.ExcelWriter(outdir+os.sep+'summary_TPM_DESeq2_DCI.xlsx')

# ==== save out the results
combined_df = pd.concat([tpm_df,deseq2_df],axis=1)
combined_df = pd.concat([combined_df,combined_dci_batch1],axis=1)
combined_df = pd.concat([combined_df,combined_dci_batch2],axis=1)

combined_df.to_csv(outdir+os.sep+'summary_TPM_DESeq2_DCI.csv')
combined_df.to_excel(writer,'All genes')

# ==== DEG
upgenes_df = combined_df[(combined_df['WT_over_Vector_padj']<0.05)&(combined_df['WT_over_Vector_log2FoldChange']>np.log2(1.25))]
dngenes_df = combined_df[(combined_df['WT_over_Vector_padj']<0.05)&(combined_df['WT_over_Vector_log2FoldChange']<-np.log2(1.25))]
# upgenes_df.shape[0]+dngenes_df.shape[0]

upgenes_df.to_csv(outdir+os.sep+'DEG_upgenes_fdr0p05_fc1p25.csv')
dngenes_df.to_csv(outdir+os.sep+'DEG_dngenes_fdr0p05_fc1p25.csv')

upgenes_df.to_excel(writer,'Up genes WToverVEC')
dngenes_df.to_excel(writer,'Down genes WToverVEC')


# ==== DCI 
k4_increased_dci_df = combined_df[(combined_df['Batch1_H3K4me3_WT_over_Vector_DCI']>5)]
k4_decreased_dci_df = combined_df[(combined_df['Batch1_H3K4me3_WT_over_Vector_DCI']<-5)]

k4_increased_dci_df.to_csv(outdir+os.sep+'Batch1_DCI_k4_increased_thre5.csv')
k4_decreased_dci_df.to_csv(outdir+os.sep+'Batch1_DCI_k4_decreased_thre5.csv')
k4_increased_dci_df.to_excel(writer,'Batch1 H3K4me3 DCI>5 WToverVEC')
k4_decreased_dci_df.to_excel(writer,'Batch1 H3K4me3 DCI<-5 WToverVEC')


# ==== DCI  batch2
k4_increased_dci_df = combined_df[(combined_df['Batch2_H3K4me3_WT_over_Vector_DCI']>5)]
k4_decreased_dci_df = combined_df[(combined_df['Batch2_H3K4me3_WT_over_Vector_DCI']<-5)]
k27_increased_dci_df = combined_df[(combined_df['Batch2_H3K27ac_WT_over_Vector_DCI']>5)]
k27_decreased_dci_df = combined_df[(combined_df['Batch2_H3K27ac_WT_over_Vector_DCI']<-5)]

k4_increased_dci_df.to_csv(outdir+os.sep+'Batch2_DCI_k4_increased_thre5.csv')
k4_decreased_dci_df.to_csv(outdir+os.sep+'Batch2_DCI_k4_decreased_thre5.csv')
k27_increased_dci_df.to_csv(outdir+os.sep+'Batch2_DCI_k27_increased_thre5.csv')
k27_decreased_dci_df.to_csv(outdir+os.sep+'Batch2_DCI_k27_decreased_thre5.csv')

k4_increased_dci_df.to_excel(writer,'Batch2 H3K4me3 DCI>5 WToverVEC')
k4_decreased_dci_df.to_excel(writer,'Batch2 H3K4me3 DCI<-5 WToverVEC')
k27_increased_dci_df.to_excel(writer,'Batch2 H3K27ac DCI>5 WToverVEC')
k27_decreased_dci_df.to_excel(writer,'Batch2 H3K27ac DCI<-5 WToverVEC')
writer.save()


