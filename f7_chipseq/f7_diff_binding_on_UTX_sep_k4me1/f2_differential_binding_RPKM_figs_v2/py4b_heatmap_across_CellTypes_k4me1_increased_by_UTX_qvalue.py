import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=10
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
                             'H3K27me3':3,
                             'MLL4':3,
                             'H3K27ac':3,
                             'H3K4me1':3,
                             'H3K4me2':3,
                             'H3K4me2APR':3,
                             'H3K4me3':3}
    
    cbar_vmax = factor_match_clim[factor]
    return cbar_vmax*0.05,cbar_vmax
        
    
              

def prepare_each_subfig(df_tmp,gs,heatmap_pos,peak_file,factor,celltype,factor_cellTypes_num):
    # read the binding pattern for each factor/celltype
    csv_file='../data_binding_patter_readCount/readCount_csv/{}_{}_on_{}_es2kb_bin200.csv'.format(celltype,factor,peak_file)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = df.loc[df_tmp.index]
    df = window_cumulative(df)
    # normalization by readcount
    norm_df = pd.read_csv('{}/f0_data_process/chip_seq/final_chipseq/sicer2_islands/get_readCount_on_sicer_islads/total_reads_in_Islands.csv'.format(project_dir),index_col=0)
    # norm_col = 'total' if re.search('UTX',factor) else 'total_in_islads'
    norm_col = 'total' 
    print(peak_file,factor,celltype,norm_col)
    norm_factor = norm_df.loc['{}_{}'.format(celltype,factor),norm_col]/1000000
    df = 50*df/norm_factor #  per kb per million mapped reads
    # return df.mean()
    pal = sns.light_palette('red',as_cmap=True)   
    vmin,vmax = return_vlim(factor,peak_file)
    # vmin=None;vmax=None
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,98))
    ax = plt.subplot(gs[0,heatmap_pos])
    g = sns.heatmap(df,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,\
                  vmin=vmin,vmax=vmax,cbar_kws={"shrink": 0.6}, rasterized=True)
    
    ax.set_ylabel('')
    cbar = g.collections[0].colorbar
    cbar.ax.set_position([.9,0.36,.8,.5])     
    cbar.set_ticks([vmin,vmax])
    cbar.set_ticklabels([0,vmax])
    cbar.ax.tick_params(labelsize=10)
    
    if heatmap_pos==0:
        ax.set_ylabel('UTX binding sites \n (n = {})\n'.format(df.shape[0]),va='baseline',fontsize=10)
    if not heatmap_pos==factor_cellTypes_num-1:
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=45,fontsize=9)
    ax.set_title('{}\n{}'.format(factor[:7], cellType_labels[celltype]),fontsize=10)
    return df.mean()





# ==== dictionary of matched colors/labels
cellType_colors = {'Vector':'tab:blue',\
                   'WT':'tab:red',\
                   'DEL':'k',\
                   'EIF':'tab:purple',\
                   'TPR':'tab:green',\
                   'MT2':'tab:purple',\
                   'FUS':'tab:gray'}

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

cellTypes = ['Vector','WT','DEL','MT2']
factors = ['UTX','H3K27me3','MLL4','H3K4me1','H3K4me2','H3K4me2APR','H3K4me3','H3K27ac']
factors = ['H3K4me1']
# peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
peak_files = ['UTX_peaks']

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
indir='../f0_data_integration_v3_RPKM/f3_combined_data_RPKM/'
# indir='../f0_data_integration_v2/f2_combined_data/'
outdir = 'f4b_heatmap_across_CellTypes_k4m1_increased_by_UTX_qvalue'
os.makedirs(outdir,exist_ok=True)

# == rank value by UTX signal 
peak_file='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed'.format(project_dir)
index_df = pd.read_csv(peak_file,sep='\t',header=None)
index_df.index = ['_'.join(i) for i in index_df.iloc[:,:3].astype(str).values]
kept_index = [i for i in index_df.index if not re.search('v',i)]
index_df = index_df.loc[kept_index]
ranked_index = index_df[8].sort_values(ascending=False).index


k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
fc_thre = 1.5
log2avg_thre = -2
avg_thres = [.25,.2,.01]

for avg_thre in avg_thres[:]:
    for peak_file in peak_files[:]:
        df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
        # ranked_index = df.UTX_WT.sort_values(ascending=False).index
        df = df.loc[ranked_index]
        # master_df = df[(df[k4me1_log2fc_col]> np.log2(fc_thre)) & (df[k4me1_log2avg_col]>log2avg_thre)] 
        master_df = df[(df[k4me1_log2fc_col]> np.log2(fc_thre)) & (df[k4me1_log2avg_col]>np.log2(avg_thre))] 
        for factor in factors[:]:
            # for each factor, check if exixt MT2 data
            factor_cellTypes_num = len([i for i in master_df.columns if i.startswith('{}_'.format(factor)) and i.endswith('FC')])+1
            factor_cellTypes = cellTypes[:factor_cellTypes_num]
            # data for composite plot
            composite_data = {}
            # heatmap of each factor
            fig = plt.figure(figsize = (.9*factor_cellTypes_num,1.3))
            width_ratio = [1]*factor_cellTypes_num
            gs = gridspec.GridSpec(1,factor_cellTypes_num,width_ratios=width_ratio,wspace=.1) 
            heatmap_pos=0
            for celltype in factor_cellTypes[:]:
                avg_binding = prepare_each_subfig(master_df,gs,heatmap_pos,peak_file,factor,celltype,factor_cellTypes_num)
                composite_data[celltype]=avg_binding
                heatmap_pos+=1
            plt.savefig(outdir+os.sep+'{}_{}_{}_binding.png'.format(avg_thre,peak_file,factor,),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
            plt.show()
            plt.close()                
            
            # == composite plot
            fig = plt.figure(figsize = (2.2,1.5))
            for celltype in factor_cellTypes[:]:
                plt.plot(composite_data[celltype], 
                         label = cellType_labels[celltype],
                         color = cellType_colors[celltype])
            plt.ylabel('RPKM'.format(factor[:7]),fontsize=11)
            plt.title('{}'.format(factor[:7]))
            # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
            plt.axes().set_yticks([0.5,1])
            plt.axes().set_xticks([0,100,200])
            plt.axes().set_xticklabels(['-2kb','0','2kb'],fontsize=10)
            plt.legend(fontsize=10,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,
                       handlelength=1,loc="upper right",
                       bbox_to_anchor=[1.01,1],
                       frameon=False)
            plt.savefig(outdir+os.sep+'composite_{}_{}_{}_binding.png'.format(avg_thre,peak_file,factor),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
            plt.show()
            plt.close()
               
    
