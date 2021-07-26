import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')



# ==== main

indir='f3b_gene_cor_diff_looping_high_threshold'
outdir = 'f6_venn_dragram'
os.makedirs(outdir,exist_ok=True)
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'

# subdirs=['f1_diff_looping','f2_diff_looping_high_thre']
factors = ['H3K4me3','H3K27ac']
loop_specificities = ['vec1_wt0','vec0_wt1']
deg_patterns = ['up','dn']
#total	
#wt_over_vec_up overlapped #wt_over_vec_up total
#wt_over_vec_dn overlapped	#wt_over_vec_dn total

for factor in factors[:]:
    csv_file='{}/{}_anchor_overlapped_genes_DEG_count.csv'.format(indir,factor)
    df = pd.read_csv(csv_file,index_col=0)
    for loop_specificity in loop_specificities[:]:
        for deg_pattern in deg_patterns[:]:
            all_genes = 20794
            total_anchor_genes = df.loc[loop_specificity,'#total']
            deg = df.loc[loop_specificity,'#wt_over_vec_{} total'.format(deg_pattern)]
            anchor_deg = df.loc[loop_specificity,'#wt_over_vec_{} overlapped'.format(deg_pattern)]
            # fisher exact test
            va = int(anchor_deg)
            vb = int(total_anchor_genes - anchor_deg)
            vc = int(deg - anchor_deg)
            vd = int(all_genes - total_anchor_genes - deg + anchor_deg)
            s,p = stats.fisher_exact([[va,vb],[vc,vd]])
            
            # ==== venn diagram
            plt.figure(figsize=(3,3))
            color1,color2 = 'tab:orange','tab:purple'
            loop_text = 'Vector' if loop_specificity=='vec1_wt0' else 'WT'
            deg_text = 'up' if deg_pattern=='up' else 'down'
            la,lb = '{}-specific loop anchor \n overlapped promoters'.format(loop_text),'{}-regulated \n genes in WT'.format(deg_text)
            out = venn2(subsets = {'10': vb, '01': vc, '11': va},set_labels=(la,lb),set_colors=(color1,color2), alpha=0.5)

            # reset the text
            label= out.get_label_by_id('11').get_text()
            out.get_label_by_id('11').set_text('{}    '.format(label))
            label= out.get_label_by_id('01').get_text()
            out.get_label_by_id('01').set_text('    {}'.format(label))
    
            # mark the text
            for text in out.set_labels:
                text.set_fontsize(16)
            for text in out.subset_labels:
                try:
                    text.set_fontsize(16)
                except:
                    pass
            plt.title(factor)
            if p<0.05:
                plt.text(x=.8,y=.75,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
            plt.savefig(outdir+os.sep+'{}_{}_{}_venn.png'.format(factor,loop_specificity,deg_pattern),bbox_inches = 'tight',pad_inches=0.1,transparent=True)
            plt.show()
            plt.close()
        

