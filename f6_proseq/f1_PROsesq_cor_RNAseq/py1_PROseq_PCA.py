import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})
# import clustering_plot
from collections import Counter
from sklearn.decomposition import PCA


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


# def get_label_color(column,label_colors_matchness):
#     for key in label_colors_matchness.keys():
#         if re.search(key,column) :
#             return [label_colors_matchness[key],key]



# 'data_1st_submit','data_20201127_1205_merged','data_20201209','data_20201220'

def get_label_marker(column):
    if re.search('_rep',column):
        return 'x'
    else:
        return 'o'
    
    
def get_label_color(column):

    label_colors_matchness={'del_cIDR':'k',\
                            'del_TPR':'purple',\
                            'Vector':'blue',\
                            'WT':'r',\
                            'UTX_FUSIDR':'g',\
                            'UTX_eIFIDR':'orange',\
                            'PCDH':'blue',\
                            'FL':'r',\
                            'DEL':'k'}
    
    for key in label_colors_matchness.keys():
        # print(key,column)
        if column.startswith(key) :
            return label_colors_matchness[key]
        
    return 'lightskyblue'
    
    
    
    
    
    
def pca_plot(df,columns,label_colors,label_markers,figname,order=None,re_labels=None,mark_text=True):

    pca=PCA(n_components=2)
    projected = pca.fit_transform(df)#;print(projected);exit()
    ratios = pca.explained_variance_ratio_
    print('pca dim:\t',projected.shape)
    
    plt.figure(figsize=(3,3))
    marker_legends,marker_list=[],[]
    color_legends=[]
    gs=[]
    for ploc in np.arange(projected.shape[0]):
        # plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        # label_mark.append(label_colors[columns[ploc]][1])
        g = plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]],marker=label_markers[columns[ploc]], s=22,label=label_markers[columns[ploc]])
        # gs.append(g)
        # plt.scatter(projected[ploc,0],projected[ploc,1],s=10,marker='^')
        # if label_colors[columns[ploc]] not in marker_list:
        #     marker_legends.append(g)
        #     marker_list.append(label_colors[columns[ploc]])
        # if mark_text:
        text = columns[ploc].split('PRO')[0]
        plt.axes().text(projected[ploc,0],projected[ploc,1],text,fontsize=11)
    # g = plt.scatter(projected[:,0],projected[:,1],c=list(label_colors.values()),marker=list(label_markers.values()), s=110)
        
    # ==== re-order the legends, if needed
    # handles,labels = plt.axes().get_legend_handles_labels()
    # if order:
        # handles,labels = re_legend_labels(handles,labels,order,re_labels)
    # plt.legend(gs)
    # handles, labels = plt.axes().legend_elements(prop="sizes", alpha=0.6)
    # plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
                # borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('PC1 ({:.1f}%)'.format(ratios[0]*100))
    plt.ylabel('PC2 ({:.1f}%)'.format(ratios[1]*100))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.show()
    plt.close()




# ==== main 
 

    
outdir = 'f1_PROseq_PCA_figs'
os.makedirs(outdir,exist_ok=True)


project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang'


## RNA-seq expression
file_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/'.format(project_dir)
tpm_file= '{}/f3_expr/tpm.csv'.format(file_dir)
with open(tpm_file) as tpm_inf:
    rna_df = pd.read_csv(tpm_inf,index_col=0)
            

## PRO-seq expression
pro_basenames= ['PCDH1PRO','PCDH2PRO','DEL31PRO','DEL32PRO','FL1PRO','FL2PRO',]
count_cols = ['RPKM','ReadCount']
dir_names = ['data_0x2_MAPQ10','data_MAPQ1']

norm_df=pd.read_excel('../data_modules_revised/Proseq_norm_factor.xlsx',sheet_name='Normalization',index_col=0)
norm_cols=['QC_dm6_q1','QC_dm6_f0x2','QC_dm6_f0x2_q10']

for dir_name in dir_names:
    pro_dir='../{}/f3_promoter_GB_UDHS_count'.format(dir_name)
    for norm_col in norm_cols:
        for count_col in count_cols:
            pro_df = pd.DataFrame()
            for pro_basename in pro_basenames:
                expr_file= '{}/{}_on_promoter.csv'.format(pro_dir,pro_basename)
                with open(expr_file) as tpm_inf:
                    df_tmp = pd.read_csv(tpm_inf,index_col=0,sep='\t')
                df_tmp = df_tmp[[count_col]].rename(columns={count_col:pro_basename})    
                # spike in normalization 
                norm_factor=norm_df.loc[pro_basename,norm_col]/1000000#;print(norm_factor)
                if count_col=='RPKM':
                    norm_factor=1
                pro_df = pd.concat([pro_df,df_tmp/norm_factor],axis=1)   

            # pro_df = pro_df[pro_df.min(axis=1)>10]
            # ==== initiate label color
            # pro_df = col_row_norm(pro_df)
            # rna_df = col_row_norm(rna_df)
            # df = pd.concat([pro_df,rna_df],axis=1).dropna()
            
            # need to be expressed in at least 5% samples
            # df = df[(df>1).sum(axis=1)>(df.shape[1]*.05)];print('df dim\t',df.shape)
            # df = np.log2(pro_df+1)
            # df = col_row_norm(pro_df)
            df = log_quantile_normalization(pro_df)
            # df = df.subtract(df.median(),axis=1)
            # df = df[[i for i in df.columns if label_colors[i]=='k']]
            # df = df[[i for i in df.columns if label_markers[i]!='x']]
            
            label_colors = {}    
            label_markers = {}
            for column in df.columns:
                label_colors[column] = get_label_color(column)   
                label_markers[column] = get_label_marker(column)   
            
                
            
            #     figname = '{}/tpm_hierarchical.pdf'.format(outdir)
            #     clustering_plot.hi_plot(np.transpose(df),df.columns,label_colors,figname)
            
            figname = '{}/QN_{}_{}_normBy_{}.pdf'.format(outdir,count_col,dir_name,norm_col)
            pca_plot(np.transpose(df),df.columns,label_colors,label_markers,figname)
            
            



