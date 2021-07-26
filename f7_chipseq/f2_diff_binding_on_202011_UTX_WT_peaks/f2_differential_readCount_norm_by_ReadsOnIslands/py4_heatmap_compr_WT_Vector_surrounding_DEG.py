import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'w'})
from matplotlib.colors import LinearSegmentedColormap
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]


    
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
        
    

def sub_heatmap_plot(df,loc,factor,celltype,deg_type,norm_col):
    
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
        # ax.set_ylabel('{} genes'.format(deg_type)) 
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=30,fontsize=14)
    ax.set_title('{}\n{}'.format(factor,celltype),fontsize=15)
              

def prepare_each_subfig(project_dir,factor,celltype,newindex,loc,genomic_region,deg_type,norm_col):
    csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_binding_pattern/rpkm_csv/{}_{}_{}_bin200_on_202011_UTX_WT_peaks.csv'.format(project_dir,celltype,factor,genomic_region)
    # csv_file='{}.csv'.format(celltype)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = window_cumulative(df)
    df = df.loc[newindex]
    # normalization by readcount
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    # norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),'total']/1000000
    # norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),'total_in_islads']/1000000
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    df = df/norm_factor
    sub_heatmap_plot(df,loc,factor,celltype,deg_type,norm_col)


outdir = 'f4_heatmap_compr_WT_Vector_surrounding_DEG'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

# rank value by UTX signal 
csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_binding_pattern/rpkm_csv/Vector_UTX_es2kb_bin200_on_202011_UTX_WT_peaks.csv'.format(project_dir)
df = pd.read_csv(csv_file,sep='\t',index_col=0)
ranked_index = df.sum(axis=1).sort_values(ascending=False).index

# heatmap compare each cell type 
factors = ['UTX','H3K27ac','H3K4me1','H3K4me2','H3K4me3','MLL4']
genomic_regions=['es2kb']
deg_types=['Up','Down','ALL',]

for norm_col in ['total','total_in_islads']:
    for deg_type in deg_types[:]:
        deg_file='{}/f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_deg_overlap/peak_overlap_promoter/peak_overlap_{}_genes_es10kb.bed'.format(project_dir,deg_type)
        deg_df = pd.read_csv(deg_file,sep='\t',index_col=3)
        newindex = [i for i in ranked_index if i in deg_df.index]
        # heatmap surrounding ALl/Up/Down genes
        for genomic_region in genomic_regions:
            for factor in factors[:]:    
                fig = plt.figure(figsize = (3,np.log10(len(newindex))))
                width_ratio = [1,1]
                gs = gridspec.GridSpec(1,2,width_ratios=width_ratio,wspace=.1) 
                # plot each cell type
                loc=0
                for celltype in ['Vector','WT']:
                    prepare_each_subfig(project_dir,factor,celltype,newindex,loc,genomic_region,deg_type,norm_col)
                    loc+=1
                        
                plt.savefig(outdir+os.sep+'{}_{}_{}_surrounding_{}genes.png'.format(norm_col,factor,genomic_region,deg_type),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
                plt.show()
                plt.close()
            
        
