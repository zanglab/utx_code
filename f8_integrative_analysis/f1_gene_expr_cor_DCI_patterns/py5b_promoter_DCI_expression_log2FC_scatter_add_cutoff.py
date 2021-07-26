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



def return_dci_df(DCI_dir,subdir,hm_mark,compr_type,suffix):

    dci_file = '{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,compr_type,suffix)
    dci_df = pd.read_csv(dci_file,sep='\t',index_col=4)
    dci_df.columns=['chr','start','end','IfOverlap','score','strand','DCI'] 
    # print(dci_df.shape)
    # dci_df.index = ['_'.join(ii) for ii in dci_df[['chr','start','end']].values.astype(str)]
    return dci_df
    

def plot_scatters(deg_df,selected_genes,compr_x,compr_y,color,label):
    x =  deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_x)].dropna()
    y =  deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_y)].dropna()
    plt.scatter(x,y,c=color,s=4,rasterized=True,label=label,alpha=1)
    return


def scatter_plot_compr_DCI(subdir,hm_mark,compr_type,suffix,dci_thre):

    compr_x =  compr_type[0]
    compr_y = compr_type[1]
    
    test_file='{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
    # print(test_file)
    if os.path.isfile(test_file):    
        dci_df_wt_over_vector = return_dci_df(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
        up_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']>dci_thre].index 
        dn_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']<-1*dci_thre].index 
        
        plt.figure(figsize=(2.7,2.7))
        plot_scatters(deg_df,deg_df.index,compr_x,compr_y,'grey','All genes')
        plot_scatters(deg_df,up_bins,compr_x,compr_y,'tab:red','Genes w/ increased DCI in WT/Vector')
        plot_scatters(deg_df,dn_bins,compr_x,compr_y,'tab:blue','Genes w/ decreased DCI in WT/Vector')
        plt.axhline(y=0,c='k',lw=1)
        plt.axvline(x=0,c='k',lw=1)
        # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
        plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,\
                   handlelength=1,loc="upper left",markerscale=3,bbox_to_anchor=[-0.1,1.32],frameon=False)
        xa,xb = cellType_labels[compr_x.split('_')[0]],cellType_labels[compr_x.split('_')[-1]]
        ya,yb = cellType_labels[compr_y.split('_')[0]],cellType_labels[compr_y.split('_')[-1]]
        plt.xlabel('Gene expression log2FC \n {} over {}'.format(xa,xb))
        plt.ylabel('Gene expression log2FC \n {} over {}'.format(ya,yb))

        plt.savefig('{}/{}/scatter_{}_{}_vs_{}{}_dci{}.png'.format(outdir,subdir,hm_mark,compr_x,compr_y,suffix,dci_thre),\
                    bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
        return up_bins,dn_bins
    return [],[]


def mark_p(p,position,val):
    if p<0.05:
        star_mark="*"
    if p<0.001:
        star_mark="**"
    if p<0.05:
        plt.text(position,val,star_mark,ha='center',va='center',fontsize=16)


def plot_box_figs(subdir,hm_mark,suffix,dci_thre,up_bins,dn_bins,num_DCI_bins_df,operator_current,operator_flag):
    

    test_file='{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
    # print(test_file)
    if os.path.isfile(test_file):    
        plt.figure(figsize=(4,3))
        position = 0  
        xticklabels = []
        for compr_col in ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']:
            # save the values in box plots
            fdr_col = '{}_padj'.format(compr_col)  
            fc_col = '{}_log2FoldChange'.format(compr_col)

            total_dci = deg_df.loc[:].dropna()
            total_dci_per = total_dci[(total_dci[fdr_col]<0.1)&(operator_current(total_dci[fc_col],float(operator_flag[1:])))]
            
            # percentage of up genes with promoter DCI pass the cutoff
            up_dci = deg_df.loc[up_bins].dropna()
            up_dci_per = up_dci[(up_dci[fdr_col]<0.1)&(operator_current(up_dci[fc_col],float(operator_flag[1:])))]
            s,p = stats.fisher_exact([[up_dci_per.shape[0],up_dci.shape[0]-up_dci_per.shape[0]],
                                      [total_dci_per.shape[0],total_dci.shape[0]-total_dci_per.shape[0]]])
            mark_p(p,position+0.4,100*up_dci_per.shape[0]/up_dci.shape[0])
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} total'.format(compr_col,'Up',operator_flag)] = up_dci.shape[0]   
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {}'.format(compr_col,'Up',operator_flag)] = up_dci_per.shape[0]     
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} s'.format(compr_col,'Up',operator_flag)] = '{:.2f}'.format(s)      
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} p'.format(compr_col,'Up',operator_flag)] = '{:.2e}'.format(p)     

 
            # percentage of down genes with promoter DCI pass the cutoff
            dn_dci = deg_df.loc[dn_bins].dropna()
            dn_dci_per = dn_dci[(dn_dci[fdr_col]<0.1)&(operator_current(dn_dci[fc_col],float(operator_flag[1:])))]
            s,p = stats.fisher_exact([[dn_dci_per.shape[0],dn_dci.shape[0]-dn_dci_per.shape[0]],
                                      [total_dci_per.shape[0],total_dci.shape[0]-total_dci_per.shape[0]]])
            mark_p(p,position+0.6,100*dn_dci_per.shape[0]/dn_dci.shape[0])
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} total'.format(compr_col,'Down',operator_flag)] = dn_dci.shape[0]   
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {}'.format(compr_col,'Down',operator_flag)] = dn_dci_per.shape[0]     
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} s'.format(compr_col,'Down',operator_flag)] = '{:.2f}'.format(s)      
            num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'{} {} {} p'.format(compr_col,'Down',operator_flag)] = '{:.2e}'.format(p)     

            

            bar_vals = [100*total_dci_per.shape[0]/total_dci.shape[0],
                        100*up_dci_per.shape[0]/up_dci.shape[0],
                        100*dn_dci_per.shape[0]/dn_dci.shape[0]]

            g=plt.bar([position+.2,position+.4,position+.6],bar_vals,width=0.15,color = ['tab:grey','tab:red','tab:blue'],linewidth=0)
            position+=1
            
            xa,xb = cellType_labels[compr_col.split('_')[0]],cellType_labels[compr_col.split('_')[-1]]    
            xticklabels.append('{} over {}'.format(xa,xb))
               

        legends = ['All genes','Genes w/ increased DCI in WT/Vector','Genes w/ decreased DCI in WT/Vector']
        plt.legend((g[0],g[1],g[2]),legends,fontsize=13,loc='upper left', bbox_to_anchor=[-0.0,1.3],\
                   frameon=False,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,handlelength=1)            
        plt.ylabel('Percentage of {} genes'.format('up' if operator_flag.startswith('>') else 'down',operator_flag),fontsize=15)    
        plt.axes().set_xticks([0.5,1.5,2.5])
        plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
        plt.axes().set_xticklabels(xticklabels,rotation=30,ha = 'right',fontsize=14)
        plt.savefig('{}/{}/box_{}_DCI{}_log2FC{}.png'.format(outdir,subdir,hm_mark,dci_thre,operator_flag),\
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

    
outdir = 'f5b_promoter_DCI_expression_log2FC_scatter_add_cutoff'
os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f2_DEG_promoter_DCI_non_normalized/f1_promoter_DCI_rename'.format(project_dir)
DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f1_DEG_promoter_DCI/f1_promoter_DCI_rename'.format(project_dir)
# DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new/bart3d_DCI_rename'.format(project_dir)
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)


subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']

compr_types = [['WT_over_Vector','DEL_over_WT'],['DEL_over_WT','EIF_over_DEL']]
hm_marks = ['H3K4me3','H3K27ac']
suffixes=['_promoter_DCI']
dci_thres = [2,5]


num_DCI_bins_df = pd.DataFrame()
for subdir in subdirs[:]: 
    outdir_tmp='{}/{}'.format(outdir,subdir)
    os.makedirs(outdir_tmp,exist_ok=True)
    for hm_mark in hm_marks[:]:
        for suffix in suffixes[:]:
            for dci_thre in dci_thres[:] :
                for compr_type in compr_types[:]:
                    up_bins,dn_bins = scatter_plot_compr_DCI(subdir,hm_mark,compr_type,suffix,dci_thre)
    
                num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'# up genes'] = len(up_bins)     
                num_DCI_bins_df.loc['{}_{}_dci{}'.format(subdir,hm_mark,dci_thre),'# dn genes'] = len(dn_bins)     
    
            
                #### box plot
                plot_box_figs(subdir,hm_mark,suffix,dci_thre,up_bins,dn_bins,num_DCI_bins_df,operator.gt,'>.1375') # np.log2(1.1)
                plot_box_figs(subdir,hm_mark,suffix,dci_thre,up_bins,dn_bins,num_DCI_bins_df,operator.lt,'<-.1375')

                # selected_bins = up_bins
                # color = 'tab:red'
                # title = 'Genes w/ increased DCI \n in WT over Vector'
                # # plot_box_figs(deg_df,selected_genes,color,title,p_thre,fc_thre,deg_count)
                # plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df)
                
                # selected_bins = dn_bins
                # color = 'tab:blue'
                # title = 'Genes w/ decreased DCI \n in WT over Vector'
                # plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df)
        

num_DCI_bins_df.to_csv(outdir+os.sep+'num_DCI_promoter_summary.csv')

    

