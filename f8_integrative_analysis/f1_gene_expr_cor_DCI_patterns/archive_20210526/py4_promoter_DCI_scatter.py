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



# def return_dci_df(DCI_dir,subdir,hm_mark,compr_type,suffix):

#     dci_file = '{}/{}/{}_{}{}.bed'.format(DCI_dir,subdir,hm_mark,compr_type,suffix)
#     dci_df = pd.read_csv(dci_file,sep='\t',header=None)
#     dci_df.columns=['chr','start','end','DCI']
#     dci_df.index = ['_'.join(ii) for ii in dci_df[['chr','start','end']].values.astype(str)]
#     return dci_df

def return_dci_df(DCI_dir,subdir,hm_mark,compr_type,suffix):

    dci_file = '{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,compr_type,suffix)
    dci_df = pd.read_csv(dci_file,sep='\t',index_col=4)
    dci_df.columns=['chr','start','end','IfOverlap','score','strand','DCI'] 
    # print(dci_df.shape)
    # dci_df.index = ['_'.join(ii) for ii in dci_df[['chr','start','end']].values.astype(str)]
    return dci_df
    

def scatter_plot_compr_DCI(subdir,hm_mark,compr_type,suffix,dci_thre):

    compr_x =  compr_type[0]
    compr_y = compr_type[1]
    
    test_file='{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,compr_x,suffix)
    # print(test_file)
    if os.path.isfile(test_file):    
        dci_df_wt_over_vector = return_dci_df(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
        up_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']>dci_thre].index 
        dn_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']<-1*dci_thre].index 
        
        dci_df_x = return_dci_df(DCI_dir,subdir,hm_mark,compr_x,suffix)
        dci_df_y = return_dci_df(DCI_dir,subdir,hm_mark,compr_y,suffix)
        
        # scatter plot
        plt.figure(figsize=(2.7,2.7))
        plt.scatter(dci_df_x.loc[:,'DCI'],dci_df_y.loc[:,'DCI'],c='tab:grey',s=3,alpha=1,rasterized=True,label='All genes')
        plt.scatter(dci_df_x.loc[up_bins,'DCI'],dci_df_y.loc[up_bins,'DCI'],c='tab:red',s=3,alpha=1,rasterized=True,label='Genes w/ increased DCI in WT/Vector')
        plt.scatter(dci_df_x.loc[dn_bins,'DCI'],dci_df_y.loc[dn_bins,'DCI'],c='tab:blue',s=3,alpha=1,rasterized=True,label='Genes w/ decreased DCI in WT/Vector')
        
        plt.axhline(y=0,c='k',lw=1)
        plt.axvline(x=0,c='k',lw=1)
        # # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
        plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,\
                    handlelength=1,loc="upper left",markerscale=3,bbox_to_anchor=[-0.1,1.32],frameon=False)
        xa,xb = cellType_labels[compr_x.split('_')[0]],cellType_labels[compr_x.split('_')[-1]]
        ya,yb = cellType_labels[compr_y.split('_')[0]],cellType_labels[compr_y.split('_')[-1]]
        plt.xlabel('Promoter DCI of {} \n {} over {}'.format(hm_mark,xa,xb))
        plt.ylabel('Promoter DCI of {} \n {} over {}'.format(hm_mark,ya,yb))
        plt.savefig('{}/{}/scatter_{}_{}_vs_{}{}_dci{}.png'.format(outdir,subdir,hm_mark,compr_x,compr_y,suffix,dci_thre),\
                    bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
        return up_bins,dn_bins
    return [],[]




def plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df):
    
    test_file='{}/{}/{}_{}{}.csv'.format(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
    
    if os.path.isfile(test_file):   
        box_vals = []
        xticklabels = []
        sig_vals,sig_colors = [],[]
        for compr_col in ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']:
            dci_df = return_dci_df(DCI_dir,subdir,hm_mark,compr_col,suffix)
            box_val = dci_df.loc[selected_bins]['DCI'].values
            # save the values in box plots
            dci_df.loc[selected_bins].to_csv('{}/{}/box_{}_{}_genes{}_dci{}_{}.csv'.format(outdir,subdir,hm_mark,title.split()[2],suffix,dci_thre,compr_col))
            s,p = stats.ttest_1samp(box_val,0)
            sig_vals.append('*' if p<0.05 else '')
            sig_colors.append('b' if s<0 else 'r')
            box_vals.append(box_val)
            xa,xb = cellType_labels[compr_col.split('_')[0]],cellType_labels[compr_col.split('_')[-1]]    
            xticklabels.append('{} over {}'.format(xa,xb))
            num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'{} {} s'.format(title.split()[2],compr_col)] = '{:.2f}'.format(s)      
            num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'{} {} p'.format(title.split()[2],compr_col)] = '{:.2e}'.format(p)     
    
            
        positions = [1,2,3]
        fig = plt.figure(figsize=(2,3))
        g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                    boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                    medianprops=dict(color='k'),showfliers=False)    
        # g = plt.violinplot(box_vals)
    
        # for position_id in np.arange(len(positions)):
        #     scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
        #     plt.scatter(scatter_x,box_vals[position_id],color=color,s=5,zorder=0,alpha=0.6,rasterized=True)
        
        # for compr_pos in [[0,1,'t'],[1,2,'t'],[2,3,'t']]:
            # mark_pvalue(compr_pos,positions,box_vals)
        plt.axes().set_xticklabels(xticklabels,rotation=30,ha='right')
        plt.ylabel('Promoter DCI of {}'.format(hm_mark),fontsize=15)
        # plt.ylim([-1,2])
        for ii in positions:
            plt.scatter(ii,np.median(box_vals[ii-1]),marker=sig_vals[ii-1],color='red',s=88)
            # plt.axes().text(ii,0,sig_vals[ii-1],fontsize=28,va='top',ha='center',color='red')
        plt.axhline(y=0,c='k',lw=1)
        plt.title(title)
        # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
        plt.savefig('{}/{}/box_{}_{}_genes{}_dci{}.png'.format(outdir,subdir,hm_mark,title.split()[2],suffix,dci_thre),\
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

    
outdir = 'f4_promoter_DCI_scatter'
os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
# DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f2_DEG_promoter_DCI_non_normalized/f1_promoter_DCI_rename'.format(project_dir)
DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f1_DEG_promoter_DCI/f1_promoter_DCI_rename'.format(project_dir)
# DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new/bart3d_DCI_rename'.format(project_dir)
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
# deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)


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
            for dci_thre in dci_thres[:]:
                for compr_type in compr_types[:]:
                    up_bins,dn_bins = scatter_plot_compr_DCI(subdir,hm_mark,compr_type,suffix,dci_thre)
    
                num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'# up genes'] = len(up_bins)     
                num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'# dn genes'] = len(dn_bins)     
    
            
                ##### box plot
                selected_bins = up_bins
                color = 'tab:red'
                title = 'Genes w/ increased DCI \n in WT over Vector'
                plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df)
                
                selected_bins = dn_bins
                color = 'tab:blue'
                title = 'Genes w/ decreased DCI \n in WT over Vector'
                plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df)
        

num_DCI_bins_df.to_csv(outdir+os.sep+'num_DCI_promoter_summary.csv')

    

