import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=13
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from matplotlib.colors import LinearSegmentedColormap
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy.interpolate import interpn
import operator


    
def window_cumulative(df,half_window=7,step=1):
    smooth_df_columns = np.arange(0,len(df.columns),step)
    smooth_df = pd.DataFrame(index=df.index,columns=smooth_df_columns)#;print(smooth_df.columns)
    for col in smooth_df_columns:
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.loc[:,col] = df.iloc[:,window_left:window_right+1].sum(axis=1)    
    #print(df,smooth_df)
    return smooth_df   


def signal_centered(df):
    center_position = int(df.shape[1]/2)
    for row in df.index:
        vals = df.loc[row]
        max_index = list(vals).index(vals.max())
        # move max to center
        if max_index<center_position:
            df.loc[row] = np.append(np.zeros(center_position-max_index),vals)[:df.shape[1]]
        elif max_index>center_position:
            df.loc[row] = np.append(vals[max_index-center_position:],np.zeros(max_index-center_position))
    return df

def window_smooth(df,half_window=7):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df


def return_vlim(norm_col):
    if norm_col=='RPKM':
        cbar_vmax = 25
    elif norm_col=='rawCount':
        cbar_vmax = 15 
    return cbar_vmax*0.06,cbar_vmax
        
    


def sub_heatmap_plot(proseq_binding_df,gs,loc,project_dir,proseq_name,norm_pattern):
    
    # selected binding pattern
    df = proseq_binding_df
    pal = sns.light_palette('red',as_cmap=True)   
    vmin,vmax = return_vlim(norm_pattern)
    # vmin=None;vmax=None
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,98))
    ax = plt.subplot(gs[loc,0])
    g = sns.heatmap(df,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,\
                  vmin=vmin,vmax=vmax,cbar_kws={"shrink": 0.5})
    
    ax.set_ylabel('')
    cbar = g.collections[0].colorbar
    cbar.set_ticks([vmin,vmax])
    cbar.set_ticklabels([0,vmax])
    ax.set_ylabel('PRO-seq peaks \n (n = {})'.format(df.shape[0]))
    ax.set_title('{}'.format(cellType_labels[proseq_names_title[proseq_name]]),fontsize=14)
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[int(len(xp)/2)],xp[-1]])
    ax.set_xticklabels(['-2kb','0','2kb'],rotation=30,fontsize=13)
    # df.to_csv(outdir+os.sep+'_{}_binding_at_UTX_with_{}_{}_pattern_{}.csv'.format(factor,flag,factor,celltype))
    return df,vmax   
    

def return_rep_combined_df(project_dir,proseq_name,rep_id,prename,norm_df,norm_pattern,norm_col):
    rep_name='{}{}'.format(proseq_name,rep_id)
    proseq_binding_file = '{}/{}_PROseq_peak_es2kb_bin200_{}_by_{}.csv'.format(indir,rep_name,norm_pattern,prename)
    proseq_binding_df = pd.read_csv(proseq_binding_file,sep='\t',index_col=0)
    # == do the normalization
    norm_factor=norm_df.loc[rep_name,norm_col]
    if norm_pattern=='rawCount':
        proseq_binding_df = 1000000*proseq_binding_df/norm_factor
    return proseq_binding_df




cellType_colors = {'Vector':'tab:blue',\
                   'WT':'tab:red',\
                   'DEL':'k',\
                   'EIF':'tab:purple',\
                   'TPR':'tab:green',\
                   'MT2':'tab:orange',\
                   'FUS':'tab:gray'}

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}



indir='f1_extract_data'
outdir = 'f4_heatmap_RPO-seq_at_enhancer_combi_rep_WT_changed'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
norm_df=pd.read_excel('{}/f6_proseq/data_modules_revised/Proseq_norm_factor.xlsx'.format(project_dir),sheet_name='Normalization',index_col=0)

proseq_names = ['PCDH','FL','DEL3',]
proseq_names_title ={'PCDH':'Vector',
                     'DEL3':'DEL',
                     'FL':'WT',}
norm_patterns=['RPKM']
prenames = ['H3K27ac_H3K4me1_with_UTX']
norm_col='total_dm6_f0x2'


diff_types = {'increased':[operator.gt,1],
             'decreased':[operator.lt,-1]}


differential_PROseq_df = pd.read_csv('f3b_RPO-seq_at_enhancer_replicate_combined_log2FC_scatter/differential_PROseq_RPKM.csv',index_col=0)
ranked_index = differential_PROseq_df.PROseq_RPKM_WT.sort_values(ascending=False).index
differential_PROseq_df = differential_PROseq_df.loc[ranked_index]

for diff_type in diff_types.keys():
    for diff_thre in [1,1.2][:1]:
        for norm_pattern in norm_patterns[:]:
            for prename in prenames[:]:   
                composite_data = {}
                for proseq_name in proseq_names[:]:
                    # == get the PROseq binding pattern from two replicates
                    proseq_binding_df1 = return_rep_combined_df(project_dir,proseq_name,'1PRO',prename,norm_df,norm_pattern,norm_col)
                    proseq_binding_df2 = return_rep_combined_df(project_dir,proseq_name,'2PRO',prename,norm_df,norm_pattern,norm_col)
                    proseq_binding_df = (proseq_binding_df1+proseq_binding_df2)/2
                    
                    compr_col = differential_PROseq_df.PROseq_RPKM_WT_over_Vector_log2FC
                    added_sign = diff_types[diff_type][1]
                    kept_df = differential_PROseq_df[diff_types[diff_type][0](compr_col,added_sign*np.log2(diff_thre))]
                    proseq_binding_df = proseq_binding_df.loc[kept_df.index]
                    # # == heatmap  
                    fig = plt.figure(figsize = (2,2.7))
                    gs = gridspec.GridSpec(1,1)#,height_ratios=height_ratios,hspace=.2) 
                    # plot each cell type
                    df,vmax = sub_heatmap_plot(proseq_binding_df,gs,0,project_dir,proseq_name,norm_pattern)
                    plt.savefig(outdir+os.sep+'WT_{}_{}_heatmap_{}_binding_by_{}_FC{}.png'.format(diff_type,norm_pattern,proseq_names_title[proseq_name],prename,diff_thre),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
                    plt.show()
                    plt.close()
                    
                    composite_data[proseq_name]=proseq_binding_df.mean()
        
                
                # == composite plot
                fig = plt.figure(figsize = (2.7,2))
                for proseq_name in proseq_names[:]:
                    celltype = proseq_names_title[proseq_name]
                    g = plt.plot(composite_data[proseq_name], 
                          label = '{}'.format(cellType_labels[celltype]),
                          color = cellType_colors[celltype],alpha=1)
                plt.ylabel('RPKM')
                plt.title('PRO-seq')
                # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
                # xp = g.get_xticks()
                # plt.axes().set_xticks([xp[0],xp[int(len(xp)/2)],xp[-1]])
                plt.axes().set_xticks([0,100,200])
                plt.axes().set_xticklabels(['-2kb','0','2kb'])
                plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.2,handletextpad=0.2,
                        handlelength=1,loc="upper right",
                        bbox_to_anchor=[1.01,1],
                        frameon=False)
                # plt.title('{}'.format(peak_file))
                plt.savefig(outdir+os.sep+'WT_{}_composite_PROseq_binding_by_{}_at_{}_FC{}.png'.format(diff_type,norm_pattern,prename,diff_thre),\
                            bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
                plt.show()
                plt.close()
                
            
        kept_df.to_csv(outdir+os.sep+'WT_{}_FC{}.csv'.format(diff_type,diff_thre))


# x = kept_df.PROseq_RPKM_WT_over_Vector_log2FC
# y = kept_df.PROseq_RPKM_DEL_over_WT_log2FC
# plt.scatter(x,y)



