import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
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
                             'H3K4me3':4}
        
    else:
        factor_match_clim = {'UTX':3,
                             'UTXFEB':3,
                             'H3K27me3':5,
                             'MLL4':5,
                             'H3K27ac':5,
                             'H3K4me1':5,
                             'H3K4me2':5,
                             'H3K4me3':5}
    
    cbar_vmax = factor_match_clim[factor]
    return cbar_vmax*0.05,cbar_vmax
        
    
 
def return_vlim_proseq(norm_col):
    if norm_col=='total_hg38':
        cbar_vmax = 1
    else:
        cbar_vmax = 40

    return cbar_vmax*0.05,cbar_vmax


             

def prepare_each_subfig(df_tmp,gs,heatmap_pos,peak_file,strand,celltype,rep,norm_col):
    # read the binding pattern for each factor/celltype
    csv_file='{}/PROseq_{}_{}_{}_on_{}_es2kb_bin200.csv'.format(proseq_pattern_dir,celltype,rep,strand,peak_file)
    df = pd.read_csv(csv_file,sep='\t',index_col=0)
    df = df.loc[df_tmp['UTX_Vector'].sort_values(ascending=False).index[:]]
    df = window_cumulative(df)
    # normalization by readcount

    norm_df=pd.read_excel('{}/f6_proseq/data_modules_revised/Proseq_norm_factor.xlsx'.format(project_dir),sheet_name='Normalization',index_col=1)
    norm_factor=norm_df.loc['{}_{}'.format(celltype,rep),norm_col]/1000000
    df = 50*df/norm_factor


    pal = sns.light_palette('red',as_cmap=True)   
    vmin,vmax = return_vlim_proseq(norm_col)
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
    if not heatmap_pos==5:
        cbar.remove()
    xp = g.get_xticks()
    ax.set_xticks([xp[0],xp[-1]])
    ax.set_xticklabels(['-2kb','2kb'],rotation=45,fontsize=13)
    strand_flag='+' if strand=='plus' else '-'
    ax.set_title('{}\n{} ({})'.format(cellType_labels[celltype],rep,strand_flag),fontsize=14)
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


project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
master_file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/f0_data_integration/f2_combined_data'.format(project_dir)
proseq_pattern_dir='../data_binding_pattern/readCount_csv/'
outdir = 'f2_heatmap_sep_replicate_with_k4me1_increased'
os.makedirs(outdir,exist_ok=True)

cellTypes = ['Vector','WT','DEL']
peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
strands = ['plus','minus']
norm_cols=['total_hg38','total_dm6_f0x2_q10','hg38_div_dm6_f0x2_q10']

k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
fc_thre = 1.5
log2avg_thre = 0


for peak_file in peak_files[:]:
    master_df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(master_file_dir,peak_file),index_col=0)
    master_df = master_df[(master_df[k4me1_log2fc_col]> np.log2(fc_thre)) & (master_df[k4me1_log2avg_col]>log2avg_thre)] 
    for strand in strands[:]:
        for norm_col in norm_cols[:2]:
            # data for composite plot
            composite_data = {}
            # heatmap of each factor
            fig = plt.figure(figsize = (6,2))
            width_ratio = [1,1,1,1,1,1]
            gs = gridspec.GridSpec(1,6,width_ratios=width_ratio,wspace=.1) 
            heatmap_pos=0
            for celltype in cellTypes[:]:
                for rep in ['rep1','rep2']:
                    avg_binding = prepare_each_subfig(master_df,gs,heatmap_pos,peak_file,strand,celltype,rep,norm_col)
                    composite_data['{}_{}'.format(celltype,rep)]=avg_binding
                    heatmap_pos+=1
            plt.savefig(outdir+os.sep+'{}_PROseq_binding_by_{}_{}.png'.format(peak_file,norm_col,strand),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
            plt.show()
            plt.close()                
            
            # == composite plot
            fig = plt.figure(figsize = (2.5,2))
            for celltype in cellTypes[:]:
                for rep in ['rep1','rep2']:
                    alpha=0.5 if rep=='rep2' else 1
                    plt.plot(composite_data['{}_{}'.format(celltype,rep)], 
                         label = '{} {}'.format(cellType_labels[celltype],rep),
                         color = cellType_colors[celltype],alpha=alpha)
            plt.ylabel('RPOseq signal \n {} strand'.format(strand))
            # plt.ylim(ymax=9 if norm_pattern=='RPKM' else 5)
            plt.axes().set_xticks([0,100,200])
            plt.axes().set_xticklabels(['-2kb','0','2kb'])
            plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.2,handletextpad=0.2,
                       handlelength=1,loc="upper right",
                       bbox_to_anchor=[1.65,1],
                       frameon=False)
            plt.title('{}'.format(peak_file))
            plt.savefig(outdir+os.sep+'composite_{}_PROseq_binding_by_{}_{}.png'.format(peak_file,norm_col,strand),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
            plt.show()
            plt.close()
           
    
