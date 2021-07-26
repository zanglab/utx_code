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
from scipy.interpolate import interpn




def return_deg(deg_df,compr_col,p_thre,fc_thre):
    log2fc_col = '{}_log2FoldChange'.format(compr_col)
    pvalue_col = '{}_padj'.format(compr_col)
    upgenes = deg_df[(deg_df[pvalue_col]<=p_thre) & (deg_df[log2fc_col]>=np.log2(fc_thre))].index
    dngenes = deg_df[(deg_df[pvalue_col]<=p_thre) & (deg_df[log2fc_col]<=-1*np.log2(fc_thre))].index
    return upgenes,dngenes
    


def plot_scatters(deg_df,selected_genes,compr_x,compr_y,color,label):
    x =  deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_x)].dropna()
    y =  deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_y)].dropna()
    plt.scatter(x,y,c=color,s=4,rasterized=True,label=label)
    return


   
########################

def plot_box_figs(deg_df,selected_genes,color,title,p_thre,fc_thre,deg_count):
    
    box_vals = []
    xticklabels = []
    sig_vals,sig_colors = [],[]
    for compr_col in ['WT_over_Vector','DEL_over_WT','EIF_over_DEL']:
        # save the values in box plots
        box_val = deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_col)].values
        deg_df.loc[selected_genes]['{}_log2FoldChange'.format(compr_col)].to_csv(outdir+os.sep+'box_{}_genes_p{}_fc{}_{}.csv'.format(title.split()[0],p_thre,fc_thre,compr_col))
        s,p = stats.ttest_1samp(box_val,0)
        sig_vals.append('*' if p<0.05 else '')
        sig_colors.append('b' if s<0 else 'r')
        box_vals.append(box_val)
        xa,xb = cellType_labels[compr_col.split('_')[0]],cellType_labels[compr_col.split('_')[-1]]    
        xticklabels.append('{} over {}'.format(xa,xb))
        deg_count.loc['p{}_fc{}'.format(p_thre,fc_thre),'{} {} s'.format(title.split()[0],compr_col)] = '{:.2f}'.format(s)
        deg_count.loc['p{}_fc{}'.format(p_thre,fc_thre),'{} {} p'.format(title.split()[0],compr_col)] = '{:.2e}'.format(p)  
        
    positions = [1,2,3]
    fig = plt.figure(figsize=(2,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    # g = plt.violinplot(box_vals)

    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=color,s=5,zorder=0,alpha=0.6,rasterized=True)
    
    # for compr_pos in [[0,1,'t'],[1,2,'t'],[2,3,'t']]:
        # mark_pvalue(compr_pos,positions,box_vals)
    plt.axes().set_xticklabels(xticklabels,rotation=30,ha='right')
    plt.ylabel('Gene expression log2FC',fontsize=15)
    # plt.ylim([-1,2])
    plt.axhline(y=0,c='k',lw=1)
    for ii in positions:
        plt.scatter(ii,np.median(box_vals[ii-1]),marker=sig_vals[ii-1],color='red',s=88)
    plt.title(title)
    # plt.legend(fontsize=16,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.savefig(outdir+os.sep+'box_{}_genes_p{}_fc{}.png'.format(title.split()[0],p_thre,fc_thre),\
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

    
project_dir="/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
project_dir="/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang"
expr_dir='{}/f0_data_process/rna_seq/data_1st_submit_STAR_RSEM_new/f6_deg/fz_deseq2_out_combined'.format(project_dir)
deg_df = pd.read_csv('{}/deseq2_combined.csv'.format(expr_dir),index_col=0)

outdir = 'f1_DEG_log2FC_scatter'
os.makedirs(outdir,exist_ok=True)

# DEG between Vector and WT
deg_count = pd.DataFrame()
for p_thre in [0.05,.01][:]:
    for fc_thre in [1.25,1.5][:]:
        upgenes,dngenes = return_deg(deg_df,'WT_over_Vector',p_thre,fc_thre)
        
        ##### scatter plot
        compr_types = [['WT_over_Vector','DEL_over_WT'],['DEL_over_WT','EIF_over_DEL']]
        for compr_type in compr_types[:]:
            compr_x =  compr_type[0]
            compr_y = compr_type[1]
            # scatter plot
            plt.figure(figsize=(2.6,2.6))
            plot_scatters(deg_df,deg_df.index,compr_x,compr_y,'grey','All genes')
            plot_scatters(deg_df,upgenes,compr_x,compr_y,'tab:red','Up genes in WT over Vector')
            plot_scatters(deg_df,dngenes,compr_x,compr_y,'tab:blue','Down genes in WT over Vector')
            
            # save and plot the correlation
            x =  deg_df['{}_log2FoldChange'.format(compr_x)].dropna()
            y =  deg_df['{}_log2FoldChange'.format(compr_y)].dropna()
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
            # output_prename = '{}_{}_{}_dci{}'.format(subdir,hm_mark,suffix,dci_thre)
            # num_DCI_bins_df.loc[output_prename,'scatter_pearsonr_s'] = r_value
            # num_DCI_bins_df.loc[output_prename,'scatter_pearsonr_p'] = p_value
            x_sort = np.sort(x)[:-1]
            plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.8)
            plt.text(.97,.97,'$r={:.2f}$ '.format(r_value),fontsize=11,transform=plt.axes().transAxes,ha='right',va='top')

            plt.axhline(y=0,c='k',lw=1)
            plt.axvline(x=0,c='k',lw=1)
            # plt.title('{} over {}'.format(cellType_labels[treatment],cellType_labels[control]))
            plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,\
                       handlelength=1,loc="upper left",markerscale=3,bbox_to_anchor=[-0.1,1.32],frameon=False)
            xa,xb = cellType_labels[compr_x.split('_')[0]],cellType_labels[compr_x.split('_')[-1]]
            ya,yb = cellType_labels[compr_y.split('_')[0]],cellType_labels[compr_y.split('_')[-1]]
            plt.xlabel('Gene expression log2FC \n {} over {}'.format(xa,xb))
            plt.ylabel('Gene expression log2FC \n {} over {}'.format(ya,yb))
            plt.savefig('{}/scatter_{}_vs_{}_p{}_fc{}.png'.format(outdir,compr_x,compr_y,p_thre,fc_thre),\
                        bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
        

        deg_count.loc['p{}_fc{}'.format(p_thre,fc_thre),'upgenes'] = len(upgenes)
        deg_count.loc['p{}_fc{}'.format(p_thre,fc_thre),'dngenes'] = len(dngenes)
        
        
        ##### box plot
        selected_genes = upgenes
        color = 'tab:red'
        title = 'Up genes in \n WT over Vector'
        plot_box_figs(deg_df,selected_genes,color,title,p_thre,fc_thre,deg_count)
        
        selected_genes = dngenes
        color = 'tab:blue'
        title = 'Down genes in \n WT over Vector'
        plot_box_figs(deg_df,selected_genes,color,title,p_thre,fc_thre,deg_count)
        
deg_count.to_csv(outdir+os.sep+'deg_count.csv')
