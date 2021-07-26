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
        cbar_vmax = 12 
    return cbar_vmax*0.06,cbar_vmax
        
    


def sub_heatmap_plot(proseq_binding_df,gs,loc,project_dir,proseq_name,utx_binding_pattern,norm_pattern):
    # get the proseq-peak ids
    peak_overlapped_dir='{}/f6_proseq/data_union_peaks/overlapped'.format(project_dir)
    proseq_id_file = '{}//union_PROseq_peak_overlapping_H3K4me1_H3K27ac_with_{}.bed'.format(peak_overlapped_dir,utx_binding_pattern)
    proseq_id_df = pd.read_csv(proseq_id_file,sep='\t',index_col=3,header=None)
    
    # selected binding pattern
    df = proseq_binding_df.loc[proseq_id_df.index]
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
    if loc==0:
        cbar.remove()
        ax.set_ylabel('PRO-seq peaks w/ \nUTX binding (#{})'.format(df.shape[0]))
        ax.set_title('{}'.format(proseq_names_title[proseq_name]),fontsize=14)
        ax.set_xticks([])
    else:
        ax.set_ylabel('PRO-seq peaks w/o \nUTX binding(#{})'.format(df.shape[0]))
        xp = g.get_xticks()
        ax.set_xticks([xp[0],xp[int(len(xp)/2)],xp[-1]])
        ax.set_xticklabels(['-2kb','0','2kb'],rotation=30,fontsize=13)
    # df.to_csv(outdir+os.sep+'_{}_binding_at_UTX_with_{}_{}_pattern_{}.csv'.format(factor,flag,factor,celltype))
    return df,vmax   
    

def return_rep_combined_df(project_dir,proseq_name,rep_id,norm_df,norm_pattern):
    rep_name='{}{}'.format(proseq_name,rep_id)
    proseq_binding_file = '{}/f6_proseq/data_0x2_MAPQ10/f4_pattern_at_PROseq_union_peak/{}_PROseq_peak_es2kb_bin200_{}.csv'.format(project_dir,rep_name,norm_pattern)
    proseq_binding_df = pd.read_csv(proseq_binding_file,sep='\t',index_col=0,header=None)
    # == do the normalization
    norm_factor=norm_df.loc[rep_name,'total_dm6_f0x2']
    if norm_pattern=='rawCount':
        proseq_binding_df = 1000000*proseq_binding_df/norm_factor
    return proseq_binding_df


outdir = 'f2_heatmap_RPO-seq_by_UTX_binding_combi_rep'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

norm_df=pd.read_excel('{}/f6_proseq/data_modules_revised/Proseq_norm_factor.xlsx'.format(project_dir),sheet_name='Normalization',index_col=0)
# proseq_names = ['PCDH1PRO','PCDH2PRO','DEL31PRO','DEL32PRO','FL1PRO','FL2PRO',]
proseq_names = ['PCDH','DEL3','FL',]
proseq_names_title ={'PCDH':'Vector',
                     'DEL3':'DEL',
                     'FL':'WT',}
utx_binding_patterns = ['bound','unbound']
norm_patterns=['RPKM','rawCount']


for proseq_name in proseq_names[:]:
    for norm_pattern in norm_patterns[:]:
        # == get the PROseq binding pattern from two replicates
        proseq_binding_df1 = return_rep_combined_df(project_dir,proseq_name,'1PRO',norm_df,norm_pattern)
        proseq_binding_df2 = return_rep_combined_df(project_dir,proseq_name,'2PRO',norm_df,norm_pattern)
        proseq_binding_df = (proseq_binding_df1+proseq_binding_df2)/2

        # == heatmap  
        fig = plt.figure(figsize = (3.5,5))
        height_ratios = [1,2]
        gs = gridspec.GridSpec(2,1,height_ratios=height_ratios,hspace=.2) 
        # plot each cell type
        df1,vmax = sub_heatmap_plot(proseq_binding_df,gs,0,project_dir,proseq_name,'WT_UTX_bound',norm_pattern)
        df2,vmax = sub_heatmap_plot(proseq_binding_df,gs,1,project_dir,proseq_name,'ALL_UTX_unbound',norm_pattern)
        plt.savefig(outdir+os.sep+'{}_heatmap_{}_binding_by.png'.format(norm_pattern,'_'.join(proseq_names_title[proseq_name].split(' '))),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.show()
        plt.close()
        
        
        # composite plot
        fig = plt.figure(figsize = (2.6,2.2))
        plt.plot(df1.mean(),label='w/ UTX',color='tab:red')
        plt.plot(df2.mean(),label='w/o UTX',color='tab:grey')
        plt.ylabel('PROseq level')
        plt.ylim(ymax=15 if norm_pattern=='RPKM' else 10)
        plt.axes().set_xticks([0,100,200])
        plt.axes().set_xticklabels(['-2kb','0','2kb'])
        plt.title(proseq_names_title[proseq_name])
        plt.legend(fontsize=12,borderaxespad=0.,labelspacing=.1,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
        plt.savefig(outdir+os.sep+'{}_composite_{}_binding_by.png'.format(norm_pattern,'_'.join(proseq_names_title[proseq_name].split(' '))),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.show()
        plt.close()
       
