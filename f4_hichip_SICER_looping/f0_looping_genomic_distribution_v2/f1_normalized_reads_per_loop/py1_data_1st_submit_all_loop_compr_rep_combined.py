import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
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

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}


# qc from MAPS
norm_file = '/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/hichip/data_1st_submit_sicer2_merged_islands_new/qc_summary_rename.xlsx'
norm_df = pd.read_excel(norm_file,index_col=1)
norm_cols = ['number_of_pairs_after_duplicate_removal',
       'number_of_intrachromosomal_pairs',
       'number_of_short-range_intrachromosomal_pairs',
       'number_of_short-range_vip_pairs',
       'number_of_long-range_intrachromosomal_pairs']


indir = '../../f0_looping_data_v2/data_1st_submission_reindex'
outdir = 'f1_data_1st_submission_rep_combined'
os.makedirs(outdir,exist_ok=True)

cellTypes = ['Vector','WT','DEL','EIF']
factors = ['H3K4me3']     
       
# prepare the box values
positions = [1,2,3,4]
for factor in factors[:]:
    for norm_col in norm_cols[:1]:
        for kept_col in ['count','expected'][:1]:            
            # plot the figs
            box_vals,xticklabels=[],[]
            for celltype in cellTypes:
                tmp_values=[] # combine the normalized values in both reps
                for rep in ['rep1','rep2']:
                    infile = '{}/{}_{}_{}.bedpe'.format(indir,factor,celltype,rep)
                    with open(infile) as inf:
                        df = pd.read_csv(inf,sep='\t',index_col=0)
                    norm_factor = norm_df.loc[norm_col,'{}_{}_{}'.format(factor,celltype,rep)]/100000000
                    normalized_vals = np.log10(df[kept_col].values/norm_factor)
                    tmp_values = np.append(tmp_values,normalized_vals)
                box_vals.append(tmp_values)
                xticklabels.append('{}'.format(cellType_labels[celltype]))
        
            # ==== plot figs
            fig = plt.figure(figsize=(3,3))
            # g = plt.violinplot(box_vals)
            g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                        boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                        medianprops=dict(color='grey'),showfliers=False)    
            
            for compr_pos in [[0,1,'t'],[1,2,'t'],[2,3,'t']]:
                mark_pvalue(compr_pos,positions,box_vals)
            plt.axes().set_xticklabels(xticklabels,rotation=45,ha='right',fontsize=15)
            plt.ylabel('log$_{{10}}$ (reads {} in loop)'.format(kept_col),fontsize=15)
            # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
            plt.title('By {}\n {}'.format(norm_col, factor, ),fontsize=15,ha='center')
            plt.savefig(outdir+os.sep+'all_loops_{}_{}_NormBy_{}.pdf'.format(factor,kept_col,norm_col),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
        
        

