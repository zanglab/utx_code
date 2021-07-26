import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats

    
def return_deg(deg_file,adjp=0.05,logfc=np.log2(1.25)):
    # read deg file, up genes
    deseq_out = pd.read_csv(deg_file)
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_upgenes = set(deseq_out_upgenes['GeneID'])
    # down genes
    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_dngenes = set(deseq_out_dngenes['GeneID'])  
    return deseq_out_upgenes,deseq_out_dngenes



def return_overlapped_gene_significance(df_count,df_names,df_index,genes,target_genes,col_name):
    overlapped = genes.intersection(target_genes)
    df_count.loc[df_index,'#{} overlapped'.format(col_name)] = len(overlapped)
    df_count.loc[df_index,'%{} overlapped'.format(col_name)] = '{:.4f}'.format(len(overlapped)/len(genes))
    df_count.loc[df_index,'#{} total'.format(col_name)] = len(target_genes)
    # fisher exact test
    s,p = stats.fisher_exact([[len(overlapped),len(genes)-len(overlapped)],[len(target_genes),20794-len(target_genes)]])
    df_count.loc[df_index,'{} odds-ratio'.format(col_name)] = s
    df_count.loc[df_index,'{} p-value'.format(col_name)] = p
    # gene names
    df_names.loc[df_index,'{}'.format(col_name)] = ','.join(overlapped)
    
    return df_count,df_names




def return_deg_overlap_info(df_count,df_names,df_index,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn):
    genes = set(i for i in genes.split(',') if len(i)>0)
    df_count.loc[df_index,'#total'] = len(genes)
    # wt_over_vec_up
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,wt_over_vec_up,'wt_over_vec_up')
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,wt_over_vec_dn,'wt_over_vec_dn')
    
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_up,'del_over_wt_up')
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_dn,'del_over_wt_dn')
    
    wt_deg = wt_over_vec_up.union(wt_over_vec_dn)
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_up.intersection(wt_deg),'del_over_wt_up_WTDEG')
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_dn.intersection(wt_deg),'del_over_wt_dn_WTDEG')
    
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_up.intersection(wt_over_vec_dn),'del_over_wt_up_WT_dn')
    df_count,df_names = return_overlapped_gene_significance(df_count,df_names,df_index,genes,del_over_wt_dn.intersection(wt_over_vec_up),'del_over_wt_dn_WT_up')
    
    return df_count,df_names






# ==== main()

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# data of the overlapped genes of each anchor for each of the union loops
indir = project_dir+os.sep+'f4_hichip_SICER_looping/f1_diff_looping/f1_diff_looping/'

outdir = 'f5_gene_cor_diff_looping_one_anchor_on_enhancer'
os.makedirs(outdir,exist_ok=True)

# get DEG comparing vector, WT, del
file_name='treated_WT_vs_ctrl_Vector.deseq2.csv'
deg_file='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out/{}'.format(project_dir,file_name)
wt_over_vec_up,wt_over_vec_dn = return_deg(deg_file)

file_name='treated_del_cIDR_vs_ctrl_WT.deseq2.csv'
deg_file='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out/{}'.format(project_dir,file_name)
del_over_wt_up,del_over_wt_dn = return_deg(deg_file)



for hm in ['H3K27ac','H3K4me3']:
    loop_file=indir+os.sep+'{}.csv'.format(hm)
    df = pd.read_csv(loop_file)
    # remove singleton loops called in only one cell type
    df = df[(df['cell_type_count']>1)&(df['cell_type_count']<5)]
    # ==== only keep those loops with both anchors on promoters
    # promoter_cols = ["anchor1_promoter_Overlapped_genes","anchor2_promoter_Overlapped_genes"]
#     df = df.dropna()
    df.loc[df.index.difference(df.dropna().index)]
    df = df.fillna('')
    # remove NAN
    df['promoter_genes'] = df["anchor1_promoter_Overlapped_genes"]+','+ df["anchor2_promoter_Overlapped_genes"]
    # for each type of diff loop, get all anchor overlapped promoters
    df = df.groupby(['cell_type'])['promoter_genes'].apply(','.join) 
    df.to_csv(outdir+os.sep+'{}_anchor_overlapped_genes.csv'.format(hm))
    
    # collect the number of genes and DEG
    df_count,df_names = pd.DataFrame(),pd.DataFrame()
    for ii in df.index:
        genes = df[ii]
        df_count,df_names = return_deg_overlap_info(df_count,df_names,ii,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn)
    
    # cell-type specific combination 
    vec_not_wt_cols = [i for i in df.index if re.search('VEC',i) and not re.search('WT',i)]
    flag = 'vec1_wt0'
    genes = ','.join([df[i] for i in vec_not_wt_cols])
    df_count,df_names = return_deg_overlap_info(df_count,df_names,flag,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn)
    
    # cell-type specific combination 
    wt_not_vec_cols = [i for i in df.index if re.search('WT',i) and not re.search('VEC',i)]
    flag = 'vec0_wt1'
    genes = ','.join([df[i] for i in wt_not_vec_cols])
    df_count,df_names = return_deg_overlap_info(df_count,df_names,flag,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn)

    # cell-type specific combination 
    wt_not_del_cols = [i for i in df.index if re.search('WT',i) and not re.search('DEL',i)]
    flag = 'wt1_del0'
    genes = ','.join([df[i] for i in wt_not_del_cols])
    df_count,df_names = return_deg_overlap_info(df_count,df_names,flag,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn)

    # cell-type specific combination 
    del_not_wt_cols = [i for i in df.index if re.search('DEL',i) and not re.search('WT',i)]
    flag = 'wt0_del1'
    genes = ','.join([df[i] for i in del_not_wt_cols])
    df_count,df_names = return_deg_overlap_info(df_count,df_names,flag,genes,wt_over_vec_up,wt_over_vec_dn,del_over_wt_up,del_over_wt_dn)
    
    df_count.to_csv(outdir+os.sep+'{}_anchor_overlapped_genes_DEG_count.csv'.format(hm))   
    df_names.to_csv(outdir+os.sep+'{}_anchor_overlapped_genes_DEG_names.csv'.format(hm))   
        
