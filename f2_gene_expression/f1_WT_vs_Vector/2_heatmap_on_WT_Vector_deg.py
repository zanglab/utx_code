import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})
# import clustering_plot
# from collections import Counter
# from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
# from sklearn.decomposition import PCA
# from sklearn.manifold import TSNE




def log_quantile_normalization(df):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = np.log2(df+1)
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized
             
def col_row_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    df = df.subtract(df.mean(axis=1),axis=0)
    return df

def log_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    return df



def return_diff_genes_deseq2(csv_file,adjp,logfc):

    deseq_out = pd.read_csv(csv_file,index_col=0)
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    # deseq_out_upgenes = set(deseq_out_upgenes['GeneID'])

    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    # deseq_out_dngenes = set(deseq_out_dngenes['GeneID'])  
    return deseq_out_upgenes,deseq_out_dngenes


# ==== main ====
if 1:    
    outdir = 'f2_heatmap'
    os.makedirs(outdir,exist_ok=True)

    file_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
    file_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'
    deseq2_file = '{}/f6_deg/f1_deseq2_out/treated_WT_vs_ctrl_Vector.deseq2.csv'.format(file_dir)
#     deseq2_file = 'treated_WT_vs_ctrl_Vector.deseq2.csv'
    adjp,logfc = 0.05,np.log2(1.25)
    upgenes,dngenes = return_diff_genes_deseq2(deseq2_file,adjp,logfc)
    deg = pd.concat([upgenes,dngenes]).sort_values(by=['log2FoldChange'],ascending=False)
    
    # ==== read TPM file
    tpm_file= '{}/f3_expr/tpm.csv'.format(file_dir)
#     tpm_file= 'tpm.csv'
    with open(tpm_file) as tpm_inf:
        df_ori = pd.read_csv(tpm_inf,index_col=0)
    print('df_ori dim\t',df_ori.shape)
    df_ori = df_ori.loc[deg.index].dropna()
    print('df_ori dim\t',df_ori.shape)#;exit()
    
    # ==== re-sort the matrix
    cols = ['Vector','WT','del_cIDR','UTX_eIFIDR','del_TPR','UTX_FUSIDR']
    xticklabels = ['Vector','WT','$\Delta$cIDR','UTX-eIF$_{IDR}$','$\Delta$TPR','UTX-FUS$_{IDR}$']
#     cols = ['WT','Vector','del_cIDR','del_TPR','UTX_eIFIDR','UTX_FUSIDR']
#     xticklabels = ['WT','Vector','$\Delta$cIDR','$\Delta$TPR','UTX_eIF$_{IDR}$','UTX_FUS$_{IDR}$']
    df = pd.DataFrame()
    for col in cols:
        # col_ids = [i for i in df_ori.columns if re.search(col,i)]
        # df[col] = df_ori[col_ids].mean(axis=1)
        col_ids = ['{}_rep1'.format(col),'{}_rep2'.format(col)]
        df = pd.concat([df,df_ori[col_ids]],axis=1)
        print(col_ids)

    df_saved = df.copy()
    # ==== normalize the matrix
    df = df[df.mean(axis=1)>0];print('df dim\t',df.shape)
    df = col_row_norm(df)
    # df = np.log2(df+1)
    
    vmax,vmin=.7,-.7
    plt.figure(figsize = (4,4))
    g = sns.heatmap(df,cmap=plt.cm.bwr,cbar=True,vmin=vmin,vmax=vmax\
                ,xticklabels=True,yticklabels=False,\
                cbar_kws={"orientation": "vertical","use_gridspec":False,'shrink':0.5}\
                )    
    
    g.collections[0].colorbar.set_label("Normalized TPM")
    g.set_ylabel('\n{} genes'.format(df.shape[0]),fontsize=22)
    xticklabels = xticklabels
    g.set_xticks([1,3,5,7,9,11])
    g.set_xticklabels(xticklabels,rotation=90,fontsize=20)
    g.vlines(np.arange(0,12,2), *g.get_ylim(),lw=1,color='w')
    g.vlines(np.arange(1,12,2), *g.get_ylim(),lw=1,color='w')
    g.xaxis.tick_top()
#     plt.show()
    plt.savefig(outdir+os.sep+'heatmap_on_DEG_WT_vs_Vector.pdf',bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()


    # ==== color bar
    fig,ax = plt.subplots(figsize=(.3,.9))
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = matplotlib.colorbar.ColorbarBase(ax,cmap=plt.cm.bwr,norm=norm,orientation='vertical')
    cb.set_ticks([vmin,0,vmax])
    ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
    plt.show()
    plt.savefig(outdir+os.sep+'colorbar.pdf',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()




#++++ save the data
    new_cols = ['{}_TPM'.format(ii) for ii in df_saved.columns]
    df_saved.columns = new_cols
    df_saved.to_csv(outdir+os.sep+'WT_over_Vector_DEG_TPM.csv')
    
    writer = pd.ExcelWriter(outdir+os.sep+'WT_over_Vector_DEG_TPM.xlsx')
    df_saved.to_excel(writer,'sheet1')
    
    #++++ add MT2 data
    mt2_deseq2_file = '/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_202011RNASEQ/star_rsem_process_human/f6_deg/f1_deseq2_out/treated_UTXWT_THP1_vs_ctrl_PCDH_THP1.deseq2.csv'
    adjp,logfc = 0.05,np.log2(1.25)
    upgenes,dngenes = return_diff_genes_deseq2(mt2_deseq2_file,adjp,logfc)
    mt2_deg = pd.concat([upgenes,dngenes]).sort_values(by=['log2FoldChange'],ascending=False)

    mt2_tpm_file='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_202011RNASEQ/star_rsem_process_human/f3_expr/tpm.csv'
    mt2_df_ori = pd.read_csv(mt2_tpm_file,index_col=0)    
    # ==== re-sort the matrix
    cols = ['PCDH_THP1','UTXWT_THP1','MT2_THP1']
    cols_matches = {'PCDH_THP1':'Vector','UTXWT_THP1':'UTX_WT','MT2_THP1':'MT2'}
    cols_renames = []
    mt2_df = pd.DataFrame()
    for col in cols:
        col_ids = ['{}_Rep1_hg38'.format(col),'{}_Rep2_hg38'.format(col)]
        cols_renames = np.append(cols_renames,['{}_rep1_TPM'.format(cols_matches[col]),'{}_rep2_TPM'.format(cols_matches[col])])
        mt2_df = pd.concat([mt2_df,mt2_df_ori[col_ids]],axis=1)
        print(col_ids)
    mt2_df.columns = cols_renames
    mt2_df = mt2_df.loc[mt2_deg.index]
    mt2_df.to_excel(writer,'sheet2')
    writer.save()
    
    
    
    
    
    
    
    
