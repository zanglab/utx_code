import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"




def return_deg(deseq_out,adjp=0.05,logfc=np.log2(1.5)):
    # read deg file, up genes
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_upgenes = deseq_out_upgenes.index
    # down genes
    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_dngenes = deseq_out_dngenes.index
    return deseq_out_upgenes,deseq_out_dngenes


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*1.01 ,1.05, 'k'
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



def plot_figs(box_vals,outdir_tmp,hm,compr_figname,norm_type):
    
    
    # plot the fig
    plt.figure(figsize=(2.6,2.6))
    
    positions=[1,2,3]
    colors = ['grey','r','b']
        
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
#         print(hm,compr_pos,s,p)

    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(['All bins','w/ H3K4me1','w/ H3K27ac'],rotation=30)    
    plt.axhline(y=0,c='k',lw=1)
    treatment,control = compr_figname.split('_')[0], compr_figname.split('_')[-1]
    plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
    # plt.xlabel('log2FoldChange')
    plt.ylabel('{} HiChIP DCI'.format(hm))
    plt.savefig(outdir_tmp+os.sep+'{}_{}_on_{}.png'.format(hm,compr_figname,norm_type),bbox_inches='tight',pad_inches=0.1,dpi=600)
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

    
project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"

outdir = 'f2_hm_binding_DCI_figs'
os.makedirs(outdir,exist_ok=True)

# peak_files = ['H3K27ac_peaks','H3K27ac_islands_increased',
#               'H3K4me1_peaks','H3K4me1_islands_increased',
#               'MLL4_peaks','MLL4_islands_increased']
subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']

name_match_df = pd.read_excel('hichip_name_matching.xlsx',index_col=0)
norm_types = ['differential_score','differential_score_after_coverageNormalization']


for subdir in subdirs[:]:
    outdir_tmp='{}/{}'.format(outdir,subdir)
    os.makedirs(outdir_tmp,exist_ok=True)
    for norm_type in norm_types[:]:
        for dci_file_basename in name_match_df.index[:]:
            print(subdir,norm_type,dci_file_basename)
            dci_file_k4 = 'f1_HM_binding_DCI/{}/{}_{}_H3K4me1_peaks_DCI.csv'.format(subdir,dci_file_basename,norm_type)
            dci_file_k27 = 'f1_HM_binding_DCI/{}/{}_{}_H3K27ac_peaks_DCI.csv'.format(subdir,dci_file_basename,norm_type)
            if os.path.isfile(dci_file_k4) and os.path.isfile(dci_file_k27):
                # print(subdir,dci_file)
                compr_figname = name_match_df.loc[dci_file_basename,'compr_figname']
                hm = name_match_df.loc[dci_file_basename,'hm']
                # read the genome wide DCI profile
                bart3d_dci_file='../f0_run_bart3d_new/{}/{}_{}.bed'.format(subdir,dci_file_basename,norm_type)
                bart3d_dci = pd.read_csv(bart3d_dci_file,header=None,sep='\t')
                genome_dci = bart3d_dci[3]
                print(compr_figname)
                # continue
                # DCI of HM overlapped regions
                k4_df = pd.read_csv(dci_file_k4,sep='\t',index_col=4)
                k4_dci = k4_df['info']
                k4_dci = [float(ii) for ele in k4_dci for ii in ele.split(',')]
                k27_df = pd.read_csv(dci_file_k27,sep='\t',index_col=4)
                k27_dci = k27_df['info']
                k27_dci = [float(ii) for ele in k27_dci for ii in ele.split(',')]
                box_vals = [genome_dci,k4_dci,k27_dci]
                # print(dci_file,[len(ii) for ii in box_vals])
                plot_figs(box_vals,outdir_tmp,hm,compr_figname,norm_type)
    

