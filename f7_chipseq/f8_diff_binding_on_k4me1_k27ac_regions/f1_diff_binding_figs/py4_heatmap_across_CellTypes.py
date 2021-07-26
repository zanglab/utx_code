import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
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
        smooth_df.loc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)    
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



def return_vlim(factor,peak_file):
    if re.search('islands',peak_file):
        factor_match_clim = {'UTX':2,
                             'UTXFEB':2,
                             'H3K27me3':4,
                             'MLL4':4,
                             'H3K27ac':4,
                             'H3K4me1':4,
                             'H3K4me2':4,
                             'H3K4me2APR':4,
                             'H3K4me3':4}
        
    else:
        factor_match_clim = {'UTX':3,
                             'UTXFEB':3,
                             'H3K27me3':5,
                             'MLL4':5,
                             'H3K27ac':5,
                             'H3K4me1':5,
                             'H3K4me2':5,
                             'H3K4me2APR':5,
                             'H3K4me3':5}
    
    cbar_vmax = factor_match_clim[factor]
    return cbar_vmax*0.05,cbar_vmax
        
    
              

def prepare_each_subfig(df_tmp,gs,heatmap_pos,peak_file,factor,celltype):
    # read the binding pattern for each factor/celltype
    csv_file='../data_binding_patter_readCount/readCount_csv/{}_{}_on_{}_es2kb_bin200.csv'.format(celltype,factor,peak_file)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = df.loc[df_tmp.index[:]]
    df = window_cumulative(df)
    # normalization by readcount
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    norm_col = 'total' if re.search('UTX',factor) else 'total_in_islads'
    print(peak_file,factor,celltype,norm_col)
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    df = 50*df/norm_factor #  per kb per million mapped reads

    pal = sns.light_palette('red',as_cmap=True)   
    vmin,vmax = return_vlim(factor,peak_file)
    # vmin=None;vmax=None
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,98))
    ax = plt.subplot(gs[0,heatmap_pos])
    g = sns.heatmap(df,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,\
                  vmin=vmin,vmax=vmax,cbar_kws={"shrink": 0.6})
    
    ax.set_ylabel('')
    cbar = g.collections[0].colorbar
    cbar.ax.set_position([.9,0.36,.8,.5])     
    cbar.set_ticks([vmin,vmax])
    cbar.set_ticklabels([0,vmax])
    
    if heatmap_pos==0:
        ax.set_ylabel('UTX binding sites \n  (#{})\n'.format(df.shape[0]),va='baseline')
    if not heatmap_pos==3:
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=45,fontsize=13)
    ax.set_title('{}\n{}'.format(factor, cellType_labels[celltype]),fontsize=14)
    return df.mean()





# ==== dictionary of matched colors/labels
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

cellTypes = ['Vector','WT','DEL','EIF']
factors = ['UTX','UTXFEB','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
peak_files = ['H3K27ac_peaks','H3K27ac_islands_increased',
              'H3K4me1_peaks','H3K4me1_islands_increased',
              'MLL4_peaks','MLL4_islands_increased']


project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
indir='../f0_data_integration/f2_combined_data/'
outdir = 'f4_heatmap_across_CellTypes'
os.makedirs(outdir,exist_ok=True)

# == rank value by UTX signal 
# csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_binding_pattern/rpkm_csv/Vector_UTX_es2kb_bin200_on_202011_UTX_WT_peaks.csv'.format(project_dir)
# index_df = pd.read_csv(csv_file,sep='\t',index_col=0)
# ranked_index = index_df.sum(axis=1).sort_values(ascending=False).index


for peak_file in peak_files[:]:
    master_df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
    for factor in factors[:]:
        # data for composite plot
        composite_data = {}
        # heatmap of each factor
        fig = plt.figure(figsize = (6,2))
        width_ratio = [1,1,1,1]
        gs = gridspec.GridSpec(1,4,width_ratios=width_ratio,wspace=.1) 
        heatmap_pos=0
        for celltype in cellTypes[:]:
            avg_binding = prepare_each_subfig(master_df,gs,heatmap_pos,peak_file,factor,celltype)
            composite_data[celltype]=avg_binding
            heatmap_pos+=1
        plt.savefig(outdir+os.sep+'{}_{}_binding.png'.format(peak_file,factor,),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.show()
        plt.close()                
        
        # == composite plot
        fig = plt.figure(figsize = (3,2))
        for celltype in cellTypes[:]:
            plt.plot(composite_data[celltype], 
                     label = cellType_labels[celltype],
                     color = cellType_colors[celltype])
        plt.ylabel('{} signal'.format(factor))
        # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
        plt.axes().set_xticks([0,100,200])
        plt.axes().set_xticklabels(['-2kb','0','2kb'])
        plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.2,handletextpad=0.2,
                   handlelength=1,loc="upper right",
                   bbox_to_anchor=[1.55,1],
                   frameon=False)
        plt.savefig(outdir+os.sep+'composite_{}_{}_binding.png'.format(peak_file,factor),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.show()
        plt.close()
           
    
