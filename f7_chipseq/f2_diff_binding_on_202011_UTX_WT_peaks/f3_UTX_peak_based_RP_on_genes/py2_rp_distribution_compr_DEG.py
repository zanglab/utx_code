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


def plot_dist(df,infile,label,color):
    df_tmp = pd.read_csv(infile,sep='\t',index_col=3,header=None)
    genes = df_tmp.index.intersection(df.index)
    kwargs = {'cumulative': True}
    sns.distplot(df.loc[genes],label=label,color=color,hist=False,hist_kws=kwargs, kde_kws=kwargs)
    # sns.distplot(df.loc[genes],label=label,color=color,hist=False)
    
    

outdir = 'f2_DEG_RP_distribution_figs'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
all_gene_file='{}/f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_deg_overlap/gene_promoters/hg38_4k_promoter_geneID.bed'.format(project_dir)
down_gene_file='{}/f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_deg_overlap/gene_promoters/WT_vs_ctrl_Vector_dngenes_promoter.bed'.format(project_dir)
up_gene_file='{}/f7_chipseq/f2_differential_binding_on_202011_UTX_WT_peaks/data_deg_overlap/gene_promoters/WT_vs_ctrl_Vector_upgenes_promoter.bed'.format(project_dir)

# rank value by UTX signal 
rp_file='gene_rp_from_UTX_WT_peak.txt'
df = pd.read_csv(rp_file,sep='\t',index_col=0,header=None)
df = np.log10(df+1)
# fig
fig = plt.figure(figsize = (3,3))
plot_dist(df,all_gene_file,'All genes','k')
plot_dist(df,up_gene_file,'Up genes','r')
plot_dist(df,down_gene_file,'Down genes','b')
# plot_dist(df,all_gene_file,'All genes','k')
# df_tmp = pd.read_csv(all_gene_file,sep='\t',index_col=3,header=None)
# genes = df_tmp.index.intersection(df.index)
# sns.distplot(df.loc[genes],label=label)

# plt.yscale('log')   
plt.ylabel('PDF')
plt.xlabel('log$_{{10}}$ RP')     
plt.legend()        
plt.savefig(outdir+os.sep+'RP_distribution.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
plt.show()
plt.close()

    
