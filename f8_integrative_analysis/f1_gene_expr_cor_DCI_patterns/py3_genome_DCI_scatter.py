import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=11
import seaborn as sns
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"



def return_dci_df(DCI_dir,subdir,hm_mark,compr_type,suffix):

    dci_file = '{}/{}/{}_{}{}.bed'.format(DCI_dir,subdir,hm_mark,compr_type,suffix)
    if os.path.isfile(dci_file):
        dci_df = pd.read_csv(dci_file,sep='\t',header=None)
        dci_df.columns=['chr','start','end','DCI']
        dci_df.index = ['_'.join(ii) for ii in dci_df[['chr','start','end']].values.astype(str)]
        return dci_df
    else:
        return None
    

def scatter_plot_compr_DCI(num_DCI_bins_df,subdir,hm_mark,compr_type,suffix,dci_thre):

    compr_x =  compr_type[0]
    compr_y = compr_type[1]
    
    test_file='{}/{}/{}_{}{}.bed'.format(DCI_dir,subdir,hm_mark,compr_y,suffix)
    print(test_file)
    if os.path.isfile(test_file):    
        dci_df_wt_over_vector = return_dci_df(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
        up_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']>dci_thre].index 
        dn_bins = dci_df_wt_over_vector[dci_df_wt_over_vector['DCI']<-1*dci_thre].index 
        
        dci_df_x = return_dci_df(DCI_dir,subdir,hm_mark,compr_x,suffix)
        dci_df_y = return_dci_df(DCI_dir,subdir,hm_mark,compr_y,suffix)
        
        # scatter plot
        plt.figure(figsize=(2.1,2.1))
        plt.scatter(dci_df_x.loc[:,'DCI'],dci_df_y.loc[:,'DCI'],c='tab:grey',s=3,alpha=1,rasterized=True,label='All bins')
        plt.scatter(dci_df_x.loc[up_bins,'DCI'],dci_df_y.loc[up_bins,'DCI'],c='tab:red',s=3,alpha=1,rasterized=True,label='Bins w/ DCI$>{}$ in WT/Vector'.format(dci_thre))
        plt.scatter(dci_df_x.loc[dn_bins,'DCI'],dci_df_y.loc[dn_bins,'DCI'],c='tab:blue',s=3,alpha=1,rasterized=True,label='Bins w/ DCI$<{}$ in WT/Vector'.format(-1*dci_thre))
        
        # save and plot the correlation
        x,y = dci_df_x.loc[:,'DCI'],dci_df_y.loc[:,'DCI']
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
        output_prename = '{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre)
        num_DCI_bins_df.loc[output_prename,'scatter_pearsonr_s'] = r_value
        num_DCI_bins_df.loc[output_prename,'scatter_pearsonr_p'] = p_value
        x_sort = np.sort(x)
        plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.8)
        plt.text(.97,.97,'$r={:.2f}$ '.format(r_value),fontsize=10,transform=plt.axes().transAxes,ha='right',va='top')

        plt.axhline(y=0,c='k',lw=1)
        plt.axvline(x=0,c='k',lw=1)
        # # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
        plt.legend(fontsize=10.5,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,\
                    handlelength=1,loc="upper left",markerscale=3,bbox_to_anchor=[-0.12,1.36],frameon=False)
        xa,xb = cellType_labels[compr_x.split('_')[0]],cellType_labels[compr_x.split('_')[-1]]
        ya,yb = cellType_labels[compr_y.split('_')[0]],cellType_labels[compr_y.split('_')[-1]]
        plt.xlabel('DCI score ({} over {})'.format(xa,xb),fontsize=12)
        plt.ylabel('DCI score ({} over {})'.format(ya,yb),fontsize=12)
        plt.savefig('{}/{}/scatter_{}_{}_vs_{}{}_dci{}.png'.format(outdir,subdir,hm_mark,compr_x,compr_y,suffix,dci_thre),\
                    bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
        return up_bins,dn_bins
    return [],[]




def plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df,flag):
    
    test_file='{}/{}/{}_{}{}.bed'.format(DCI_dir,subdir,hm_mark,'WT_over_Vector',suffix)
    
    if os.path.isfile(test_file):   
        box_vals = []
        xticklabels = []
        sig_vals,sig_colors = [],[]
        for compr_col in ['WT_over_Vector','DEL_over_WT','EIF_over_DEL','TPR_over_WT']:
            dci_df = return_dci_df(DCI_dir,subdir,hm_mark,compr_col,suffix)
            if dci_df is not None:
                box_val = dci_df.loc[selected_bins]['DCI'].values
                dci_df.loc[selected_bins].to_csv('{}/{}/box_{}_{}_bins{}_dci{}_{}.csv'.format(outdir,subdir,hm_mark,flag,suffix,dci_thre,compr_col))
                s,p = stats.ttest_1samp(box_val,0)
                sig_vals.append('*' if p<0.05 else '')
                sig_colors.append('b' if s<0 else 'r')
                box_vals.append(box_val)
                xa,xb = cellType_labels[compr_col.split('_')[0]],cellType_labels[compr_col.split('_')[-1]]    
                xticklabels.append('{} over {}'.format(xa,xb))
                num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'{} {} s'.format(title.split()[2],compr_col)] = '{:.2f}'.format(s)      
                num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'{} {} p'.format(title.split()[2],compr_col)] = '{:.2e}'.format(p)     
    
            
        positions = np.arange(len(box_vals))
        fig = plt.figure(figsize=(.46*len(box_vals),2.2))
        g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                    boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                    medianprops=dict(color='grey'),showfliers=False)    
        # g = plt.violinplot(box_vals)
    
        # for position_id in np.arange(len(positions)):
        #     scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
        #     plt.scatter(scatter_x,box_vals[position_id],color=color,s=5,zorder=0,alpha=0.6,rasterized=True)
        
        # for compr_pos in [[0,1,'t'],[1,2,'t'],[2,3,'t']]:
            # mark_pvalue(compr_pos,positions,box_vals)
        plt.axes().set_xticklabels(xticklabels,rotation=30,ha='right',fontsize=12)
        plt.ylabel('DCI score'.format(hm_mark),fontsize=12)
        # plt.ylim([-1,2])
        for ii in positions:
            plt.scatter(ii,np.median(box_vals[ii]),marker=sig_vals[ii],color='red',s=77)
        plt.axhline(y=0,c='k',lw=1)
        plt.title(title,fontsize=12)
        # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
        plt.savefig('{}/{}/box_{}_{}_bins{}_dci{}.png'.format(outdir,subdir,hm_mark,flag,suffix,dci_thre),\
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

    
outdir = 'f3_genome_DCI_scatter'
os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
DCI_dir='{}/f5_hichip/f1_hichip_bart3d_new/f0_run_bart3d_new/bart3d_DCI_rename'.format(project_dir)
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/f1_deseq2_out'.format(project_dir)
# expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
# deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)


subdirs=['bart3d_dis200k_data_1st_submit','bart3d_dis200k_data202008',
         'bart3d_dis500k_data_1st_submit','bart3d_dis500k_data202008']

compr_types = [['WT_over_Vector','DEL_over_WT'],['DEL_over_WT','EIF_over_DEL'],['WT_over_Vector','TPR_over_WT']]
hm_marks = ['H3K4me3','H3K27ac']
suffixes=['_differential_score','_differential_score_after_coverageNormalization']
dci_thres = [5,7,9]


num_DCI_bins_df = pd.DataFrame()
for subdir in subdirs[:]: 
    outdir_tmp='{}/{}'.format(outdir,subdir)
    os.makedirs(outdir_tmp,exist_ok=True)
    for hm_mark in hm_marks[:]:
        for suffix in suffixes[1:]:
            for dci_thre in dci_thres[:1]:
                for compr_type in compr_types[:]:
                    up_bins,dn_bins = scatter_plot_compr_DCI(num_DCI_bins_df,subdir,hm_mark,compr_type,suffix,dci_thre)
    
                    if compr_type[1]=='DEL_over_WT':
                        num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'# up bins'] = len(up_bins)     
                        num_DCI_bins_df.loc['{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre),'# dn bins'] = len(dn_bins)     
            
                    
                        ##### box plot
                        selected_bins = up_bins
                        color = 'tab:red'
                        title = 'Bins w/ DCI$>{}$ \n in WT over Vector'.format(dci_thre)
                        plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df,'increased')
                        
                        selected_bins = dn_bins
                        color = 'tab:blue'
                        title = 'Bins w/ DCI$<{}$ \n in WT over Vector'.format(-1*dci_thre)
                        plot_box_figs(subdir,hm_mark,suffix,selected_bins,color,title,dci_thre,num_DCI_bins_df,'decreased')
            

num_DCI_bins_df.to_csv(outdir+os.sep+'num_DCI_bins_summary.csv')

    

