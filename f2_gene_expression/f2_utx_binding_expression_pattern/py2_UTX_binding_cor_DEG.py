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
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),90)*1.01 ,1.05, 'k'
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



def plot_figs(box_vals,outdir,compr_type,dis,peak_file):
    
    
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


    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(['All genes','w/ UTX binding','w/ UTX binding and \n increased k4me1'],rotation=45)    
    plt.axhline(y=0,c='k',lw=1)
    treatment,control = compr_type.split('_')[0], compr_type.split('_')[-1]
    plt.title('{} \n {} over {}'.format(peak_file, cellType_labels[treatment],cellType_labels[control]))
    plt.ylabel('log2FoldChange')
    plt.savefig('{}/{}_dis{}bp_compr_{}.png'.format(outdir,peak_file,dis,compr_type,),bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.show()
    plt.close()
    
    
def utx_nearby_gene_log2FC(utx_df,target_df,deg_df):
    utx_df = utx_df.loc[target_df.index].dropna()
    overlapped_genes = utx_df[utx_df['IfOverlap']==1]['overlapped_genes']
    overlapped_genes = [gene for ele in overlapped_genes for gene in ele.split(',')]
    overlapped_genes = set(overlapped_genes).intersection(deg_df.index)
    overlapped_log2FC = deg_df.loc[overlapped_genes]['log2FoldChange'].dropna()
    return overlapped_log2FC



# ==== main() 

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}

deseq2_labels = {'WT_over_Vector':'treated_WT_vs_ctrl_Vector.deseq2.csv',\
                 'DEL_over_WT':'treated_del_cIDR_vs_ctrl_WT.deseq2.csv',\
                 'EIF_over_DEL':'treated_UTX_eIFIDR_vs_ctrl_del_cIDR.deseq2.csv'}

project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
master_file_dir='{}/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/f0_data_integration/f2_combined_data'.format(project_dir)
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)

outdir = 'f2_utx_binding_DEG_figs'
os.makedirs(outdir,exist_ok=True)

peak_files = ['UTX_peaks','UTX_islands','UTXFEB_islands','UTXFEB_peaks']
k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
fc_thres = [1.5,2]
fc_thre = 1.5
log2avg_thre = 0

extend_dis = [0,2000,10000,50000]
compr_types = ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']

# expression of genes near UTX           
for peak_file in peak_files[:]:
    # read the master dataframe
    master_df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(master_file_dir,peak_file),index_col=0)
    master_df_tmp = master_df[(master_df[k4me1_log2fc_col]> np.log2(fc_thre)) & (master_df[k4me1_log2avg_col]>log2avg_thre)] 
    for compr_type in compr_types[:]:
        deg_file = expr_dir+os.sep+deseq2_labels[compr_type]
        deg_df = pd.read_csv(deg_file,index_col=0)
        # upgenes,dngenes = return_deg(deg_df)
        # get the genes surrounding UTX
        for dis in extend_dis[:]:
            utx_file = 'f1_UTX_binding_promoter_overlap/{}_tss_es{}bp.csv'.format(peak_file,dis)
            utx_df = pd.read_csv(utx_file,sep='\t',index_col=4)
            all_log2FC =  deg_df['log2FoldChange'].dropna()
            utx_log2FC = utx_nearby_gene_log2FC(utx_df,master_df,deg_df)
            utx_k4_log2FC = utx_nearby_gene_log2FC(utx_df,master_df_tmp,deg_df)
            box_vals = [all_log2FC,utx_log2FC,utx_k4_log2FC]
            try:
                plot_figs(box_vals,outdir,compr_type,dis,peak_file)
            except:
                pass
            utx_log2FC.to_csv('{}/_{}_dis{}bp_compr_{}_log2FC.csv'.format(outdir,peak_file,dis,compr_type))
            utx_k4_log2FC.to_csv('{}/_{}_dis{}bp_compr_{}_log2FC_with_increased_k4me1.csv'.format(outdir,peak_file,dis,compr_type))
    

