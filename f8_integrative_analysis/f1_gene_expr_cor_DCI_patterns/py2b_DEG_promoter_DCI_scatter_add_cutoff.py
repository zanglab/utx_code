import sys,argparse,operator
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



def return_deg(deg_df,compr_col,p_thre,fc_thre):
    log2fc_col = '{}_log2FoldChange'.format(compr_col)
    pvalue_col = '{}_padj'.format(compr_col)
    upgenes = deg_df[(deg_df[pvalue_col]<=p_thre) & (deg_df[log2fc_col]>=np.log2(fc_thre))].index
    dngenes = deg_df[(deg_df[pvalue_col]<=p_thre) & (deg_df[log2fc_col]<=-1*np.log2(fc_thre))].index
    return upgenes,dngenes

def scatter_plot_compr_DCI(deg_df,outdir_tmp,hm_mark,compr_type,p_thre,fc_thre):

    upgenes,dngenes = return_deg(deg_df,'WT_over_Vector',p_thre,fc_thre)
    compr_x =  compr_type[0]
    compr_y = compr_type[1]
    dci_file_x='{}/{}/{}_{}_promoter_DCI.csv'.format(DCI_dir,subdir,hm_mark,compr_x)
    dci_file_y='{}/{}/{}_{}_promoter_DCI.csv'.format(DCI_dir,subdir,hm_mark,compr_y)
    if os.path.isfile(dci_file_x):
        dci_df_x = pd.read_csv(dci_file_x,sep='\t',index_col=4)
        dci_df_y = pd.read_csv(dci_file_y,sep='\t',index_col=4)
        
        # scatter plot
        plt.figure(figsize=(2.7,2.7))
        plt.scatter(dci_df_x.loc[:,'info'],dci_df_y.loc[:,'info'],c='tab:grey',s=3,rasterized=True,label='All genes')
        plt.scatter(dci_df_x.loc[upgenes,'info'],dci_df_y.loc[upgenes,'info'],c='tab:red',s=3,rasterized=True,label='Up genes in WT over Vector')
        plt.scatter(dci_df_x.loc[dngenes,'info'],dci_df_y.loc[dngenes,'info'],c='tab:blue',s=3,rasterized=True,label='Down genes in WT over Vector')
        
        plt.axhline(y=0,c='k',lw=1)
        plt.axvline(x=0,c='k',lw=1)
        # # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
        plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,\
                    handlelength=1,loc="upper left",markerscale=3,bbox_to_anchor=[-0.1,1.32],frameon=False)
        xa,xb = cellType_labels[compr_x.split('_')[0]],cellType_labels[compr_x.split('_')[-1]]
        ya,yb = cellType_labels[compr_y.split('_')[0]],cellType_labels[compr_y.split('_')[-1]]
        plt.xlabel('Promoter DCI of {} \n {} over {}'.format(hm_mark,xa,xb))
        plt.ylabel('Promoter DCI of {} \n {} over {}'.format(hm_mark,ya,yb))
        plt.savefig('{}/scatter_{}_{}_vs_{}_p{}_fc{}.png'.format(outdir_tmp,hm_mark,compr_x,compr_y,p_thre,fc_thre),\
                    bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
    
    return upgenes,dngenes




def return_dci_df(DCI_dir,subdir,hm_mark,compr_type,suffix):

    dci_file = '{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,compr_type,suffix)
    dci_df = pd.read_csv(dci_file,sep='\t',index_col=4)
    dci_df.columns=['chr','start','end','IfOverlap','score','strand','DCI'] 
    return dci_df


def mark_p(p,position,val):
    if p<0.05:
        star_mark="*"
    if p<0.001:
        star_mark="**"
    if p<0.05:
        plt.text(position,val,star_mark,ha='center',va='center',fontsize=16)


def plot_box_figs(subdir,hm_mark,suffix,upgenes,dngenes,num_DCI_bins_df,operator,operator_flag):
    
    test_file='{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
    
    if os.path.isfile(test_file):   
        plt.figure(figsize=(4,3))
        position = 0  
        xticklabels = []
        for compr_col in ['WT_over_Vector','DEL_over_WT','EIF_over_DEL'][:]:
            dci_df = return_dci_df(DCI_dir,subdir,hm_mark,compr_col,suffix)
            
            # percentage of DCI pass a cutoff
            total_dci = dci_df.loc[:].dropna()
            total_dci_per = total_dci[operator(total_dci['DCI'],int(operator_flag[1:]) )]
            
            # percentage of up genes with promoter DCI pass the cutoff
            up_dci = dci_df.loc[upgenes].dropna()
            up_dci_per = up_dci[operator(up_dci['DCI'],int(operator_flag[1:]) )]
            s,p = stats.fisher_exact([[up_dci_per.shape[0],up_dci.shape[0]-up_dci_per.shape[0]],
                                      [total_dci_per.shape[0],total_dci.shape[0]-total_dci_per.shape[0]]])
            mark_p(p,position+0.4,100*up_dci_per.shape[0]/up_dci.shape[0])
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} total'.format(compr_col,'Up',operator_flag)] = up_dci.shape[0]   
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {}'.format(compr_col,'Up',operator_flag)] = up_dci_per.shape[0]     
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} s'.format(compr_col,'Up',operator_flag)] = '{:.2f}'.format(s)      
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} p'.format(compr_col,'Up',operator_flag)] = '{:.2e}'.format(p)     

            
            # percentage of down genes with promoter DCI pass the cutoff
            dn_dci = dci_df.loc[dngenes].dropna()
            dn_dci_per = dn_dci[operator(dn_dci['DCI'],int(operator_flag[1:]) )]
            s,p = stats.fisher_exact([[dn_dci_per.shape[0],dn_dci.shape[0]-dn_dci_per.shape[0]],
                                      [total_dci_per.shape[0],total_dci.shape[0]-total_dci_per.shape[0]]])
            mark_p(p,position+0.6,100*dn_dci_per.shape[0]/dn_dci.shape[0])
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} total'.format(compr_col,'Down',operator_flag)] = dn_dci.shape[0]   
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {}'.format(compr_col,'Down',operator_flag)] = dn_dci_per.shape[0]     
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} s'.format(compr_col,'Down',operator_flag)] = '{:.2f}'.format(s)      
            num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'{} {} {} p'.format(compr_col,'Down',operator_flag)] = '{:.2e}'.format(p)     
            
            # plot the bar
            bar_vals = [100*total_dci_per.shape[0]/total_dci.shape[0],
                        100*up_dci_per.shape[0]/up_dci.shape[0],
                        100*dn_dci_per.shape[0]/dn_dci.shape[0]]

            g=plt.bar([position+.2,position+.4,position+.6],bar_vals,width=0.15,color = ['tab:grey','tab:red','tab:blue'],linewidth=0)
            position+=1
            
            xa,xb = cellType_labels[compr_col.split('_')[0]],cellType_labels[compr_col.split('_')[-1]]    
            xticklabels.append('{} over {}'.format(xa,xb))
               
        legends = ['All genes','Up genes in WT over Vector','Down genes in WT over Vector']
        plt.legend((g[0],g[1],g[2]),legends,fontsize=13,loc='upper left', bbox_to_anchor=[-0.0,1.3],\
                   frameon=False,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,handlelength=1)            
        plt.ylabel('Percentage of DCI${}$'.format(operator_flag),fontsize=15)    
        plt.axes().set_xticks([0.5,1.5,2.5])
        plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
        plt.axes().set_xticklabels(xticklabels,rotation=30,ha = 'right',fontsize=14)
        plt.savefig('{}/{}/bar_{}_p{}_fc{}_DCI{}.png'.format(outdir,subdir,hm_mark,p_thre,fc_thre,operator_flag),\
                    bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
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

    
outdir = 'f2b_DEG_promoter_DCI_scatter_add_cutoff'
os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f2_DEG_promoter_DCI_non_normalized/f1_promoter_DCI_rename'.format(project_dir)
# DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f1_DEG_promoter_DCI/f1_promoter_DCI_rename'.format(project_dir)
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)


subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']

compr_types = [['WT_over_Vector','DEL_over_WT'],['DEL_over_WT','EIF_over_DEL']]
hm_marks = ['H3K4me3','H3K27ac']

num_DCI_bins_df = pd.DataFrame()
for subdir in subdirs[:]: 
    outdir_tmp='{}/{}'.format(outdir,subdir)
    os.makedirs(outdir_tmp,exist_ok=True)
    for hm_mark in hm_marks[:]:
        for p_thre in [0.05,.01][:1]:
            for fc_thre in [1.25,1.5][:1]:
                for compr_type in compr_types[:]:
                    upgenes,dngenes = scatter_plot_compr_DCI(deg_df,outdir_tmp,hm_mark,compr_type,p_thre,fc_thre)                    

                num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'upgenes'] = len(upgenes)
                num_DCI_bins_df.loc['{}_{}_p{}_fc{}'.format(subdir,hm_mark,p_thre,fc_thre),'dngenes'] = len(dngenes)
        
                suffix = '_promoter_DCI' 
                plot_box_figs(subdir,hm_mark,suffix,upgenes,dngenes,num_DCI_bins_df,operator.gt,'>2')
                plot_box_figs(subdir,hm_mark,suffix,upgenes,dngenes,num_DCI_bins_df,operator.lt,'<-2')


num_DCI_bins_df.to_csv(outdir+os.sep+'num_DCI_promoter_summary.csv')
