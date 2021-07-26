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
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
# import clustering_plot
# from collections import Counter
# from sklearn.decomposition import PCA    
    
    

def box_plot(df,figname,xticklabels,count_col):
    
    # prepare the box values
    box_vals=[]
    positions = [1,2,3,4,5,6]
    colors=['b','b','r','r','k','k']
    for col in df.columns:
        box_vals.append(df[col].values)


    # ==== plot figs
    fig = plt.figure(figsize=(3.5,3))
    # g = plt.violinplot(box_vals)
    g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.04,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=1,zorder=0,alpha=0.7,rasterized=True)
    
    # for compr_pos in [[0,1,'t'],[0,2,'b']]:
    #     mark_pvalue(compr_pos,positions,box_vals)
    
    plt.axes().set_xticklabels(xticklabels,rotation=45,ha='right')
    plt.ylabel('log$_{{10}}$ {}'.format(count_col),fontsize=13)
    # plt.title(re_names[compr_name_ii],fontsize=13)
    # plt.axhline(y=0,c='k',lw=1)
    # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()


# ==== main 
 

    
outdir = 'f2_PROseq_on_HK_figs'
os.makedirs(outdir,exist_ok=True)


project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
            

## PRO-seq expression
pro_basenames= ['PCDH1PRO','PCDH2PRO','FL1PRO','FL2PRO','DEL31PRO','DEL32PRO',]
xticklabels= ['Vector rep1','Vector rep2','WT rep1','WT rep2','DEL rep1','DEL rep2',]
count_cols = ['ReadCount','RPKM',]
dir_names = ['data_0x2_MAPQ10']
genomic_regions=['GB','promoter']

hk_genes = pd.read_csv('../data_modules_revised/HK_genes.txt',sep=r'\s+',index_col=0,header=None)
norm_df=pd.read_excel('../data_modules_revised/Proseq_norm_factor.xlsx',sheet_name='Normalization',index_col=0)
norm_cols=['total_hg38','total_dm6_f0x2_q10','hg38_div_dm6_f0x2_q10']

for dir_name in dir_names[:1]:
    pro_dir='../{}/f3_promoter_GB_UDHS_count'.format(dir_name)
    for norm_col in norm_cols[:]:
        for count_col in count_cols[:]:
            for genomic_region in genomic_regions[:2]:
                pro_df = pd.DataFrame()
                for pro_basename in pro_basenames:
                    expr_file= '{}/{}_on_{}.csv'.format(pro_dir,pro_basename,genomic_region)
                    with open(expr_file) as tpm_inf:
                        df_tmp = pd.read_csv(tpm_inf,index_col=0,sep='\t')
                    df_tmp = df_tmp[[count_col]].rename(columns={count_col:pro_basename})   
                    df_tmp = df_tmp.loc[hk_genes.index].dropna()
                    # spike in normalization 
                    norm_factor=norm_df.loc[pro_basename,norm_col]/1000000#;print(norm_factor)
                    if count_col=='RPKM':
                        norm_factor=1
                    pro_df = pd.concat([pro_df,df_tmp/norm_factor],axis=1)   
                # min_thre = 1 if count_col=='RPKM' else 10
                # df = pro_df[pro_df.min(axis=1)>min_thre];print(df.shape)
                # df = pro_df.loc[pro_df.sum(axis=1).sort_values(ascending=False)[:200000].index]
                df = pro_df.loc[pro_df.var(axis=1).sort_values(ascending=False)[:200000].index]
                df = np.log10(df+1)
                # plot the figs
                figname = '{}/{}_{}_{}_normBy_{}.pdf'.format(outdir,count_col,genomic_region,dir_name,norm_col)
                box_plot(df,figname,xticklabels,count_col)
            
            



