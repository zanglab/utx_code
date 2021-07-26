import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"




def return_deg(deseq_out,adjp=0.05,logfc=np.log2(1.25)):
    # read deg file, up genes
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_upgenes = deseq_out_upgenes.index
    # down genes
    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_dngenes = deseq_out_dngenes.index
    return deseq_out_upgenes,deseq_out_dngenes


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),94)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),4)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    s_label='dn' if s>0 else 'up'
    p_label='{:.0e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05:
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, '{} {}'.format(p_label,s_label), ha='center', va='bottom', color=col,fontsize=11)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.65, y2*.8, '{} {}'.format(p_label,s_label), ha='center', va='top', color=col,fontsize=11)
    return s,p



def plot_figs(deg_file,dci_file,outdir_tmp,hm,compr_figname):
    
    # read info
    deg_df = pd.read_csv(deg_file,index_col=0)
    dci_df = pd.read_csv(dci_file,sep='\t',index_col=4)
    # keep only those shared genes
    shared_genes= deg_df.index.intersection(dci_df.index)
    deg_df = deg_df.loc[shared_genes]
    dci_df = dci_df.loc[shared_genes]
    # select DEG
    upgenes,dngenes = return_deg(deg_df)
    print(outdir_tmp,hm,compr_figname, '#up-genes',len(upgenes), '#down-genes', len(dngenes))
    # plot the fig
    plt.figure(figsize=(2.6,2.6))
    
    positions=[1,2,3]
    colors = ['grey','r','b']
    
    box_vals = []
    for genes in [shared_genes,upgenes,dngenes]:
        box_vals.append(dci_df.loc[genes,'info'].values)
    # plt.violinplot(box_vals)
    
    g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
                
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=True)


    for compr_pos in [[0,1,'t'],[0,2,'b']]:
        s,p = mark_pvalue(compr_pos,positions,box_vals)
        #print(hm,compr_pos,s,p)

    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(['All genes','Up genes','Down genes'],rotation=45)    
    plt.axhline(y=0,c='k',lw=1)
    treatment,control = compr_figname.split('_')[0], compr_figname.split('_')[-1]
    plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
    # plt.xlabel('log2FoldChange')
    plt.ylabel('DCI score'.format(hm))
    plt.savefig(outdir_tmp+os.sep+'{}_{}.png'.format(hm,compr_figname),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    
    

# ==== main() 

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

outdir = 'f2_expr_DCI_figs'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)

subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']


name_match_df = pd.read_excel('dci_deg_matching.xlsx',index_col=0)

for subdir in subdirs[:]:
    outdir_tmp='{}/{}'.format(outdir,subdir)
    os.makedirs(outdir_tmp,exist_ok=True)
    promoter_dci_files=glob.glob('f1_promoter_DCI/{}/*.csv'.format(subdir))
    for dci_file in promoter_dci_files:
        dci_file_basename = os.path.basename(dci_file)
        if dci_file_basename in name_match_df.index:
            # print(dci_file_basename)
            compr_figname = name_match_df.loc[dci_file_basename,'compr_figname']
            hm = name_match_df.loc[dci_file_basename,'hm']
            deg_file = expr_dir+os.sep+name_match_df.loc[dci_file_basename,'deg_file']
            plot_figs(deg_file,dci_file,outdir_tmp,hm,compr_figname)
    

