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
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*.8 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.0e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    p_label = '{} {}'.format(p_label,'up' if s<0 else 'down')
    print(s,p)
    print(p_label)
    # if p<0.05:
    #     if compr_pos[2] == 't':
    #         plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y, y*h, y*h, y], lw=1, c=col)
    #         plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=10)
    #     else:
    #         plt.plot([x1*1.05, x1*1.05, x2*0.95, x2*0.95], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
    #         plt.text((x1+x2)*.5, y2*.9, p_label, ha='center', va='top', color=col,fontsize=10)
    return s,p,p_label



## ==== main 

# re_names = ['Vector','WT','$\Delta$cIDR','UTX_eIF$_{IDR}$']
compr_names=['K4M3WT_over_K4M3PCDH','K4M3DEL3_over_K4M3WT','K4M3EIF_over_K4M3DEL3',
             'K27ACWT_over_K27ACVEC','K27ACDEL_over_K27ACWT','K27ACEIF_over_K27ACDEL']
re_names = ['H3K4me3 WT over Vector','H3K4me3 $\Delta$cIDR over WT','H3K4me3 UTX_eIF$_{IDR}$ over $\Delta$cIDR',
            'H3K27ac WT over Vector','H3K27 $\Delta$cIDR over WT','H3K27ac UTX_eIF$_{IDR}$ over $\Delta$cIDR']


compr_names=['K4M3WT_over_K4M3PCDH',
             'K27ACWT_over_K27ACVEC',]
re_names = ['H3K4me3 WT over Vector',
            'H3K27ac WT over Vector',]


for bart_dis in ['200k']:
    indir = 'overlapped_{}'.format(bart_dis)
    outdir = 'f2_figs_compr_dci_{}'.format(bart_dis)
    os.makedirs(outdir,exist_ok=True)
    for compr_name_ii in np.arange(len(compr_names))[:]:
        compr_name = compr_names[compr_name_ii]
        all_dci='../data_dci_{}/{}_differential_score.bed'.format(bart_dis,compr_name)
        increased_file='{}/{}_dci_overlap_sicer_increased.bed'.format(indir,compr_name)
        decreased_file='{}/{}_dci_overlap_sicer_decreased.bed'.format(indir,compr_name)

        # prepare the box values
        box_vals=[]
        positions = [1,2,3]
        colors=['grey','r','b']
        for infile in [all_dci,increased_file]:
            with open(infile) as inf:
                df = pd.read_csv(inf,sep='\t',header=None)
            df.columns = ['chr','start','end','score']
            box_vals.append(df['score'].values)
        unbound_value = np.setdiff1d(box_vals[0],box_vals[1])
        box_vals.append(unbound_value)
        

        # ==== plot figs
        fig = plt.figure(figsize=(2.6,2.6))
        # g = plt.violinplot(box_vals)
        g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                    boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                    medianprops=dict(color='grey'),showfliers=False)    
        
        # for position_id in np.arange(len(positions)):
        #     scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
        #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=7,zorder=0,alpha=0.99,rasterized=True)
        
        for compr_pos in [[0,1,'t'],[0,2,'b'],[1,2,'b']]:
            s,p,p_label = mark_pvalue(compr_pos,positions,box_vals)
            print(compr_pos,s,p,p_label)
        plt.axes().set_xticklabels(['All bins','UTX bound','UTX unbound'],rotation=30,ha='right')
        plt.ylabel('DCI',fontsize=13)
        plt.title(re_names[compr_name_ii],fontsize=13)
        plt.axhline(y=0,c='k',lw=1)
        # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
        plt.savefig(outdir+os.sep+'{}_{}_compr_dci.pdf'.format(bart_dis,compr_name),bbox_inches='tight',pad_inches=0.1,dpi=600)
        plt.show()
        plt.close()
    


