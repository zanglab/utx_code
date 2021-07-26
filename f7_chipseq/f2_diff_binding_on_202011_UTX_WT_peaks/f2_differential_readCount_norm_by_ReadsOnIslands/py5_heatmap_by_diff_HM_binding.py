import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
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

def return_vlim(factor,norm_col):
    if norm_col=='total':
        factor_match_clim = {'UTX':1.5,
                             'MLL4':1.5,
                             'H3K27ac':2,
                             'H3K4me1':2,
                             'H3K4me2':2,
                             'H3K4me3':2}

    elif norm_col=='total_in_islads':
        factor_match_clim = {'UTX':5,
                             'MLL4':4,
                             'H3K27ac':4,
                             'H3K4me1':4,
                             'H3K4me2':2,
                             'H3K4me3':4}
    
    cbar_vmax = factor_match_clim[factor]
    return cbar_vmax*0.06,cbar_vmax
        
    

def sub_heatmap_plot(gs,df,loc,factor,celltype,norm_col,flag):
    
    pal = sns.light_palette('red',as_cmap=True)   
    # vmin=cbarvmin,vmax = vmax,
    # vmin=.2;vmax=4
    vmin,vmax = return_vlim(factor,norm_col)
    # vmin=None;vmax=None
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,98))
    ax = plt.subplot(gs[0,loc])
    g = sns.heatmap(df,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,\
                  vmin=vmin,vmax=vmax,cbar_kws={"shrink": 0.5})
    
    ax.set_ylabel('')
    cbar = g.collections[0].colorbar
    # cbar.set_clim(vmax*.05,vmax)
    cbar.set_ticks([vmin,vmax])
    cbar.set_ticklabels([0,vmax])
    
    if loc==0:
        ax.set_ylabel('UTX binding sites w/ \n {} {} (#{})\n'.format(flag,factor,df.shape[0]),va='baseline')
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=30,fontsize=13)
    ax.set_title('{}\n{}'.format(factor,celltype),fontsize=14)
              

def prepare_each_subfig(gs,project_dir,factor,celltype,newindex,loc,norm_col,flag):
    csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_binding_pattern/rpkm_csv/{}_{}_es2kb_bin200_on_202011_UTX_WT_peaks.csv'.format(project_dir,celltype,factor)
    # csv_file='{}.csv'.format(celltype)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = window_cumulative(df)
    df = df.loc[newindex]
    # normalization by readcount
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    df = df/norm_factor
    sub_heatmap_plot(gs,df,loc,factor,celltype,norm_col,flag)
    # df.to_csv(outdir+os.sep+'_{}_binding_at_UTX_with_{}_{}_pattern_{}.csv'.format(factor,flag,factor,celltype))
    return df


def heatmap_compr_diff_binding(diff_df,flag,ranked_index,outdir,factor,treatment,control,fc_cutoff):

    # == select the index/regions for heatmap
    newindex = [i for i in ranked_index if i in diff_df.index]
    norm_col = 'total_in_islads'
    
    # == heatmap by differential HM binding
    fig = plt.figure(figsize = (3,2))
    width_ratio = [1,1]
    gs = gridspec.GridSpec(1,2,width_ratios=width_ratio,wspace=.1) 
    # plot each cell type
    # loc=0
    # for celltype in ['Vector','WT']:
    #     prepare_each_subfig(gs,project_dir,factor,celltype,newindex,loc,norm_col,flag)
    #     loc+=1
    df1 = prepare_each_subfig(gs,project_dir,factor,'Vector',newindex,0,norm_col,flag)
    df2 = prepare_each_subfig(gs,project_dir,factor,'WT',newindex,1,norm_col,flag)
                
    plt.savefig(outdir+os.sep+'fc{}_{}_binding_at_UTX_with_{}_{}.png'.format(fc_cutoff,factor,flag,factor),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()
    # diff_df.to_csv(outdir+os.sep+'_fc{}_{}_binding_at_UTX_with_{}_{}_index.csv'.format(fc_cutoff,factor,flag,factor))
            
    
    # == composite plot
    fig = plt.figure(figsize = (2.5,2))
    plt.plot(df1.mean(),label='Vector',color='tab:grey')
    plt.plot(df2.mean(),label='WT',color='tab:red')
    plt.ylabel('{} level'.format(factor))
    # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
    plt.axes().set_xticks([0,100,200])
    plt.axes().set_xticklabels(['-2kb','0','2kb'])
    plt.title('{} {}'.format(factor, flag),fontsize=13)
    plt.legend(fontsize=12,borderaxespad=0.,labelspacing=.1,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.savefig(outdir+os.sep+'composite_fc{}_{}_binding_at_UTX_with_{}_{}.png'.format(fc_cutoff,factor,flag,factor),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()
       



outdir = 'f5_heatmap_compr_WT_Vector_by_diff_HM_binding'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

# == rank value by UTX signal 
csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_binding_pattern/rpkm_csv/Vector_UTX_es2kb_bin200_on_202011_UTX_WT_peaks.csv'.format(project_dir)
index_df = pd.read_csv(csv_file,sep='\t',index_col=0)
ranked_index = index_df.sum(axis=1).sort_values(ascending=False).index

# == read the master file with all different binding info
master_csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/fz_data_combination/f2_combined_data/combined_differential_NormReadCount.csv'.format(project_dir)
master_df = pd.read_csv(master_csv_file,index_col=0)

# == heatmap compare each cell type 
# genomic_regions=['es2kb']
# norm_col = ['total','total_in_islads']
# factors = ['UTX','H3K27ac','H3K4me1','H3K4me2','H3K4me3','MLL4']
factors = ['H3K4me1','H3K27ac',]
celltype_pairs = [['WT','Vector'],]
fc_cutoffs=[1.5,2]

for fc_cutoff in fc_cutoffs:
    diff_thre_lower=-1*np.log2(fc_cutoff)
    diff_thre_higher=np.log2(fc_cutoff)
    
    for factor in factors[:]:
        for celltype_pair in celltype_pairs[:]:
            treatment = '{}_{}'.format(factor,celltype_pair[0])
            control = '{}_{}'.format(factor,celltype_pair[1])
            diff_logavg = master_df['{}_over_{}_log2AVG'.format(treatment,control)]
            diff_logfc = master_df['{}_over_{}_log2FC'.format(treatment,control)]
            # == get regions with differential HM binding
            up_df = master_df[(diff_logavg>0)&(diff_logfc>=diff_thre_higher)]
            middle_df = master_df[(diff_logavg>0)&(diff_logfc>diff_thre_lower)&(diff_logfc<diff_thre_higher)]
            down_df = master_df[(diff_logavg>0)&(diff_logfc<=diff_thre_lower)]
            heatmap_compr_diff_binding(up_df,'increased',ranked_index,outdir,factor,treatment,control,fc_cutoff)
            heatmap_compr_diff_binding(middle_df,'not-changed',ranked_index,outdir,factor,treatment,control,fc_cutoff)
            heatmap_compr_diff_binding(down_df,'decreased',ranked_index,outdir,factor,treatment,control,fc_cutoff)
            # newindex = [i for i in ranked_index if i in deg_df.index]
