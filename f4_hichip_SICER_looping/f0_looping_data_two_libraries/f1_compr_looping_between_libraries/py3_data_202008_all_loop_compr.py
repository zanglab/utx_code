import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.rcParams['font.size']=18
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*.98 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.0e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    p_label = '{}\n{}'.format(p_label,'up' if s<0 else 'down')
    if p<0.05:
        if compr_pos[2] == 't':
            plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=10)
        else:
            plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=10)




## ==== main 

indir = '../data_202008_reindex'
outdir = 'f3_data_202008_all_loop'
os.makedirs(outdir,exist_ok=True)

sample_names = ['K4M3PCDH','K4M3WT','K4M3DEL3','K4M3EIF','K4M3DTPR']
# number_of_pairs_after_duplicate_removal
norm_factors = {'K4M3PCDH': 99417566 ,'K4M3WT':97428856,'K4M3DEL3':115072544,'K4M3EIF':114846745,'K4M3DTPR':133855777}
# number_of_long-range_intrachromosomal_pairs
# norm_factors = {'K4M3PCDH': 35896339 ,'K4M3WT':38208305,'K4M3DEL3':40579723,'K4M3EIF':44971872,'K4M3DTPR':41177211}
re_names = ['Vector','WT','$\Delta$cIDR','UTX_eIF$_{IDR}$','$\Delta$TPR']


# prepare the box values
positions = [1,2,3,4,5]
for kept_col in ['count','expected']:            
    # plot the figs
    box_vals=[]
    for sample_name in sample_names:
        infile = '{}/{}.bedpe'.format(indir,sample_name)
        with open(infile) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=0)
            normalized_vals = np.log10(100000000*df[kept_col].values/norm_factors[sample_name]+1)
            box_vals.append(normalized_vals)


    # ==== plot figs
    fig = plt.figure(figsize=(4,3))
    # g = plt.violinplot(box_vals)
    g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    
#     for compr_pos in [[0,1,'t'],[1,2,'t'],[2,3,'t'],[3,4,'t']]:
#         mark_pvalue(compr_pos,positions,box_vals)
    plt.axes().set_xticklabels(re_names,rotation=30)
    plt.ylabel('log$_{{10}}$ (reads {} in loop)'.format(kept_col),fontsize=15)
    # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.savefig(outdir+os.sep+'all_loops_{}.pdf'.format(kept_col),bbox_inches='tight',pad_inches=0.1,dpi=600)
#     plt.show()
    plt.close()



