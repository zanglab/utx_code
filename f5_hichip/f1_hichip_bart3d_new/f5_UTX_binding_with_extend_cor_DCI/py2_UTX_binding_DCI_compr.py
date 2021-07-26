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



def plot_figs(box_vals,outdir_tmp,hm,compr_figname,peak_file,extend_dis):
    
    
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


    for compr_pos in [[0,1,'t'],[1,2,'t'],[0,2,'b']]:
        s,p = mark_pvalue(compr_pos,positions,box_vals)
#         print(hm,compr_pos,s,p)

    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(['All bins','w/ UTX binding','w/ UTX binding and \n increased k4me1'],rotation=30)    
    plt.axhline(y=0,c='k',lw=1)
    treatment,control = compr_figname.split('_')[0], compr_figname.split('_')[-1]
    plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
    # plt.xlabel('log2FoldChange')
    plt.ylabel('{} HiChIP DCI'.format(hm))
    plt.savefig(outdir_tmp+os.sep+'{}_{}_on_{}_es{}.png'.format(hm,compr_figname,peak_file,extend_dis),bbox_inches='tight',pad_inches=0.1,dpi=600)
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
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)

outdir = 'f2_utx_binding_DCI_figs'
os.makedirs(outdir,exist_ok=True)

master_file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/f0_data_integration/f2_combined_data'.format(project_dir)
peak_files = ['UTX_peaks','UTX_islands','UTXFEB_islands','UTXFEB_peaks']
k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
fc_thres = [1.5,2]
fc_thre = 1.5
log2avg_thre = 0


subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']
name_match_df = pd.read_excel('dci_utx_matching.xlsx',index_col=0)


for peak_file in peak_files[:]:
    # read the master dataframe
    master_df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(master_file_dir,peak_file),index_col=0)
    master_df_tmp = master_df[(master_df[k4me1_log2fc_col]> np.log2(fc_thre)) & (master_df[k4me1_log2avg_col]>log2avg_thre)] 
    # for each hichip data type and each utx binding file
    for subdir in subdirs[:]:
        outdir_tmp='{}/{}'.format(outdir,subdir)
        os.makedirs(outdir_tmp,exist_ok=True)
        for dci_file_basename in name_match_df.index[:]:
            for extend_dis in [0,20000,50000]:
                dci_file = 'f1_UTX_binding_DCI/{}/{}_{}_es{}_DCI.csv'.format(subdir,dci_file_basename,peak_file,extend_dis)
                if os.path.isfile(dci_file):
                    print(subdir,dci_file)
                    utx_df = pd.read_csv(dci_file,sep='\t',index_col=4)
                    compr_figname = name_match_df.loc[dci_file_basename,'compr_figname']
                    hm = name_match_df.loc[dci_file_basename,'hm']
                    # read the genome wide DCI profile
                    bart3d_dci_file='../f0_run_bart3d_new/{}/{}_differential_score_after_coverageNormalization.bed'.format(subdir,dci_file_basename)
                    bart3d_dci = pd.read_csv(bart3d_dci_file,header=None,sep='\t')
                    # get the dci of all bins/ utx-binding regions/ those with increased k4
                    genome_dci = bart3d_dci[3]
                    utx_binding_dci = utx_df.loc[master_df.index].dropna()['info'].values
                    utx_binding_dci = [float(ii) for ele in utx_binding_dci for ii in ele.split(',')]
                    utx_with_increased_k4 = utx_df.loc[master_df_tmp.index].dropna()['info'].values
                    utx_with_increased_k4 = [float(ii) for ele in utx_with_increased_k4 for ii in ele.split(',')]
                    
                    box_vals = [genome_dci,utx_binding_dci,utx_with_increased_k4]
                    print(dci_file,[len(ii) for ii in box_vals])
                    plot_figs(box_vals,outdir_tmp,hm,compr_figname,peak_file,extend_dis)
    

