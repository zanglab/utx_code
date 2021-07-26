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
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from matplotlib.colors import LinearSegmentedColormap
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy.interpolate import interpn
from scipy import stats



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*.7 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.0e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    p_label = '{}\n{}'.format(p_label,'up' if s<0 else 'down')
    
    if p<0.05:
        if compr_pos[2] == 't':
            plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=13)
        else:
            plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.9, p_label, ha='center', va='top', color=col,fontsize=13)






indir = 'f5_heatmap_compr_WT_Vector_by_diff_HM_binding'
outdir ='f6_boxes_compr_WT_Vector_by_diff_HM_binding'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'


# == heatmap compare each cell type 
# genomic_regions=['es2kb']
# norm_col = ['total','total_in_islads']
# factors = ['UTX','H3K27ac','H3K4me1','H3K4me2','H3K4me3','MLL4']
factors = ['H3K4me1','H3K27ac',]
celltype_pairs = [['WT','Vector'],]
diff_thre_lower=-1*np.log2(1.5)
diff_thre_higher=np.log2(1.5)
diff_patterns = ['increased','decreased','not-changed']

# == read the master file with all different binding info
master_csv_file='{}//f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/fz_data_combination/f2_combined_data/combined_differential_NormReadCount.csv'.format(project_dir)
master_df = pd.read_csv(master_csv_file,index_col=0)


for factor in factors[:]:
    for diff_pattern in diff_patterns[:3]:
        # get the index of UTX binding sites
        index_file='{}/_{}_binding_at_UTX_with_{}_{}_index.csv'.format(indir,factor,diff_pattern,factor)
        index_df = pd.read_csv(index_file,index_col=0)
        # get the HM binding signals from the master df
        box_vals = []
        master_col1='{}_{}'.format(factor,'Vector')
        master_col2='{}_{}'.format(factor,'WT')
        vals1 = master_df.loc[index_df.index][master_col1]
        vals2 = master_df.loc[index_df.index][master_col2]
        box_vals.append(np.log2(vals1+.1))
        box_vals.append(np.log2(vals2+.1))

        positions = [1,2]
        colors=['grey','grey']

        # ==== plot figs
        fig = plt.figure(figsize=(2.6,2.6))
        # g = plt.violinplot(box_vals)
        g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                    boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                    medianprops=dict(color='grey'),showfliers=False)    
        
        for position_id in np.arange(len(positions)):
            scatter_x = np.random.normal(positions[position_id],0.05,len(box_vals[position_id]))
            plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=6,zorder=0,alpha=0.99,rasterized=True)
        
        for compr_pos in [[0,1,'t']]:
            mark_pvalue(compr_pos,positions,box_vals)
        plt.axes().set_xticklabels(['Vector','WT'],rotation=30,ha='right')
        plt.ylabel('{} signal'.format(factor),fontsize=16)
        # plt.title(re_names[compr_name_ii],fontsize=13)
        # plt.axhline(y=0,c='k',lw=1)
        # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
        plt.savefig(outdir+os.sep+'{}_binding_at_UTX_with_{}_{}.pdf'.format(factor,diff_pattern,factor),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()


