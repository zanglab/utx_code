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
        ax.set_ylabel('{} {} (#{})\n'.format(flag,factor,df.shape[0]),va='baseline')
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=30,fontsize=13)
    ax.set_title('{}\n{}'.format(factor,celltype),fontsize=14)
              

def prepare_each_subfig(gs,project_dir,factor,celltype,newindex,loc,norm_col,flag):
    csv_file='{}/f7_chipseq/f3_differential_binding_on_sicer_df//data_binding_pattern/readCount_csv/{}_{}_es2kb_bin200.csv'.format(project_dir,celltype,factor)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = df.loc[newindex]
    df = window_cumulative(df)
    # normalization by readcount
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    df = df/norm_factor
    sub_heatmap_plot(gs,df,loc,factor,celltype,norm_col,flag)
    # df.to_csv(outdir+os.sep+'_{}_binding_at_UTX_with_{}_{}_pattern_{}.csv'.format(factor,flag,factor,celltype))
    return df


def heatmap_compr_diff_binding(diff_df,flag,outdir,factor,norm_col,fc_cutoff):

    newindex = diff_df.index
    prename='{}_fc{}_{}_binding_WT_vs_Vector_{}'.format(norm_col,fc_cutoff,factor,flag)
    # == heatmap by differential HM binding
    fig = plt.figure(figsize = (3,2))
    width_ratio = [1,1]
    gs = gridspec.GridSpec(1,2,width_ratios=width_ratio,wspace=.1) 
    df1 = prepare_each_subfig(gs,project_dir,factor,'Vector',newindex,0,norm_col,flag)
    df2 = prepare_each_subfig(gs,project_dir,factor,'WT',newindex,1,norm_col,flag)
    plt.savefig(outdir+os.sep+'{}_heatmap.png'.format(prename),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()
    diff_df.to_csv(outdir+os.sep+'_{}_index.csv'.format(prename))
            
    
    # == composite plot
    fig = plt.figure(figsize = (2.5,2.2))
    plt.plot(df1.mean(),label='Vector',color='tab:grey')
    plt.plot(df2.mean(),label='WT',color='tab:red')
    plt.ylabel('{} level'.format(factor))
    # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
    plt.axes().set_xticks([0,100,200])
    plt.axes().set_xticklabels(['-2kb','0','2kb'])
    plt.title('{} {}'.format(factor, flag),fontsize=13)
    plt.legend(fontsize=12,borderaxespad=0.,labelspacing=.1,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.savefig(outdir+os.sep+'{}_composite.png'.format(prename),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()
       



outdir = 'f1_heatmap_compr_WT_Vector'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

# == heatmap compare each cell type 
# factors = ['UTX','H3K27ac','H3K4me1','H3K4me2','H3K4me3','MLL4']
factors = ['H3K4me1','H3K27ac',]
norm_cols = ['total_in_islads','total']
fc_cutoffs = [1.5,2]

for factor in factors[:]:
    for norm_col in norm_cols[:]:
        for fc_cutoff in fc_cutoffs[:]:
            sicer_df_results='{}/f7_chipseq/f3_differential_binding_on_sicer_df/data_sicer_df/sicer_df_results/WT_vs_Vector_{}/WT_{}_treat-and-Vector_{}_treat-W200-G600-summary'.format(project_dir,factor,factor,factor)
            sicer_df = pd.read_csv(sicer_df_results,sep='\t')
            # the id in binding pattern starts from 1
            sicer_df.index = sicer_df.index+1
            up_df = sicer_df[sicer_df['Fc_A_vs_B']>=fc_cutoff]
            down_df = sicer_df[sicer_df['Fc_B_vs_A']>=fc_cutoff]
            middle_df = sicer_df[(sicer_df['Fc_A_vs_B']<fc_cutoff)&(sicer_df['Fc_B_vs_A']<fc_cutoff)]
            # print(up_df.shape[0]+down_df.shape[0]+middle_df.shape[0])    
            heatmap_compr_diff_binding(up_df,'increased',outdir,factor,norm_col,fc_cutoff)
            heatmap_compr_diff_binding(down_df,'decreased',outdir,factor,norm_col,fc_cutoff)
            heatmap_compr_diff_binding(middle_df,'not-changed',outdir,factor,norm_col,fc_cutoff)





