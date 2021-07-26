import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=13
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from matplotlib.colors import LinearSegmentedColormap
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy.interpolate import interpn
from scipy import stats


    

def compr_diff_proseq_scatter(df,stats_df,outdir,factor,count_type,compr_type,batch,hm):
    

    treatment_x = compr_pairs[ii][0]
    control_x = compr_pairs[ii][1]
    # compr_pairs[ii+1]
    treatment_y = compr_pairs[ii+1][0]
    control_y = compr_pairs[ii+1][1]
    factor_logFC_x = '{}_{}_over_{}_{}_log2FC'.format(factor,treatment_x,factor,control_x)
    factor_logFC_y = '{}_{}_over_{}_{}_log2FC'.format(factor,treatment_y,factor,control_y)
    
    x = df[factor_logFC_x]
    y = df[factor_logFC_y]
            

    treatment,control = compr_type.split('_over_')
    dci_col = '{}_{}_{}_DCI'.format(batch,hm,compr_type)
    rpkm_col = '{}_{}_{}_log2FC'.format(count_type,factor,compr_type)
    
    if dci_col in df.columns:
        x = df[dci_col]
        y = df[rpkm_col]
        
        # == stats test
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
        output_prename = '{}_{}_cor_DCI_{}_{}_{}'.format(factor,count_type,compr_type,batch,hm)
        stats_df.loc[output_prename,'pearsonr_s'] = r_value
        stats_df.loc[output_prename,'pearsonr_p'] = p_value
    
    
        plt.figure(figsize=(3,3))
        data, x_e,y_e = np.histogram2d(x,y,bins=20,density = True)
        z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
        z[np.where(np.isnan(z))] = 0.0
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
        # plt.scatter(x,y,rasterized=True,s=7,c='k')
        plt.axhline(y=0,c='gray',lw=.8,ls='--')
        plt.axvline(x=0,c='gray',lw=.8,ls='--')
        plt.xlabel('{} {} $\Delta$DCI ({}/{})'.format(batch,hm,treatment,control))
        plt.ylabel('{} $\Delta${} ({}/{})'.format(factor,count_type, treatment,control))
        # plt.title('$\Delta${}'.format(factor))
        x_sort = np.sort(x)
        plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.9)
        # plt.text(.13,1.04,'$r={:.2f}, p={:.1e}$'.format(r_value,p_value),fontsize=13,transform=plt.axes().transAxes)
        plt.text(.97,.97,'$r={:.2f}$'.format(r_value),fontsize=13,transform=plt.axes().transAxes,ha='right',va='top')
        plt.savefig(outdir+os.sep+'{}.pdf'.format(output_prename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
        plt.show()
        plt.close()

    

def return_differential_score(df,treatment,control,outdir):
    # == log2FC
    pc=0.1
    t = 'PROseq_RPKM_{}'.format(treatment)
    c = 'PROseq_RPKM_{}'.format(control)
    x = np.log2(df[[t,c]].mean(axis=1)+pc)
    y = np.log2((df[t]+pc)/(df[c]+pc))
    df['PROseq_RPKM_{}_over_{}_log2AVG'.format(treatment,control)] = x
    df['PROseq_RPKM_{}_over_{}_log2FC'.format(treatment,control)] = y
    
    
 
    
def return_rep_combined_df(project_dir,proseq_name,rep_id,prename,norm_pattern):
    rep_name='{}{}'.format(proseq_name,rep_id)
    proseq_binding_file = '{}/{}_PROseq_peak_es2kb_bin200_{}_by_{}.csv'.format(indir,rep_name,norm_pattern,prename)
    proseq_binding_df = pd.read_csv(proseq_binding_file,sep='\t',index_col=0)
    # == do the normalization
    # norm_factor=norm_df.loc[rep_name,norm_col]
    # if norm_pattern=='rawCount':
        # proseq_binding_df = 1000000*proseq_binding_df/norm_factor
    return proseq_binding_df




cellType_colors = {'Vector':'tab:blue',\
                   'WT':'tab:red',\
                   'DEL':'k',\
                   'EIF':'tab:purple',\
                   'TPR':'tab:green',\
                   'MT2':'tab:orange',\
                   'FUS':'tab:gray'}

cellType_labels= {'Vector':'Vector',\
                  'WT':'WT',\
                  'DEL':'$\Delta$cIDR',\
                  'EIF':'UTX-eIF$_{IDR}$',\
                  'TPR':'$\Delta$TPR',\
                  'MT2':'MT2',\
                  'FUS':'UTX-FUS$_{IDR}$'}



indir='f1_extract_data'
outdir = 'f3_RPO-seq_at_enhancer_replicate_combined_log2FC_scatter'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
# norm_df=pd.read_excel('{}/f6_proseq/data_modules_revised/Proseq_norm_factor.xlsx'.format(project_dir),sheet_name='Normalization',index_col=0)
# norm_col='total_dm6_f0x2'

proseq_names = ['PCDH','FL','DEL3',]
proseq_names_title ={'PCDH':'Vector',
                     'DEL3':'DEL',
                     'FL':'WT',}
# norm_patterns=['RPKM','rawCount']
# prenames = ['H3K4me1_overlapping_H3K27ac_overlapping_UTX','H3K27ac_H3K4me1_with_UTX']
norm_patterns=['RPKM']
prenames = ['H3K27ac_H3K4me1_with_UTX']

proseq_RPKM_df = pd.DataFrame()
for norm_pattern in norm_patterns[:]:
    for prename in prenames[:]:   
        composite_data = {}
        for proseq_name in proseq_names[:]:
            # == get the PROseq binding pattern from two replicates
            proseq_binding_df1 = return_rep_combined_df(project_dir,proseq_name,'1PRO',prename,norm_pattern)
            proseq_binding_df2 = return_rep_combined_df(project_dir,proseq_name,'2PRO',prename,norm_pattern)
            proseq_binding_df = (proseq_binding_df1+proseq_binding_df2)/2
            proseq_RPKM_df['PROseq_RPKM_{}'.format(proseq_names_title[proseq_name])] = proseq_binding_df.mean(axis=1)

# get the differential PROseq signals 
return_differential_score(proseq_RPKM_df,'WT','Vector',outdir)            
return_differential_score(proseq_RPKM_df,'DEL','WT',outdir)            
# return_differential_score(proseq_RPKM_df,'EIF','DEL',outdir)            
proseq_RPKM_df.to_csv(outdir+os.sep+'differential_PROseq_RPKM.csv')

compr_pairs = [['WT','Vector'],['DEL','WT']]
stats_df = pd.DataFrame()
for ii in np.arange(len(compr_pairs)-1)[:]:
            # compr_pairs[ii]
            treatment_x = compr_pairs[ii][0]
            control_x = compr_pairs[ii][1]
            # compr_pairs[ii+1]
            treatment_y = compr_pairs[ii+1][0]
            control_y = compr_pairs[ii+1][1]
            logFC_x = 'PROseq_RPKM_{}_over_{}_log2FC'.format(treatment_x,control_x)
            logFC_y = 'PROseq_RPKM_{}_over_{}_log2FC'.format(treatment_y,control_y)
            
            x = proseq_RPKM_df[logFC_x]
            y = proseq_RPKM_df[logFC_y]
            
            # == stats test
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)    
            s,p = stats.pearsonr(x, y)    
            output_prename = 'PROseq_RPKM_{}_over_{}_vs_{}_over_{}'.format(treatment_x,control_x,treatment_y,control_y)
            stats_df.loc[output_prename,'pearsonr_s'] = s
            stats_df.loc[output_prename,'pearsonr_p'] = p
            stats_df.loc[output_prename,'r_value'] = r_value
            stats_df.loc[output_prename,'p_value'] = p_value

            plt.figure(figsize=(2.5,2.5))
            data, x_e,y_e = np.histogram2d(x,y,bins=30,density = True)
            z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
            z[np.where(np.isnan(z))] = 0.0
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=6,rasterized=True)
            # plt.scatter(x,y,rasterized=True,s=7,c='k')
            plt.axhline(y=0,c='gray',lw=.8,ls='--')
            plt.axvline(x=0,c='gray',lw=.8,ls='--')
            plt.xlabel('$\Delta$PROseq ({} over {})'.format( cellType_labels[treatment_x],cellType_labels[control_x]))
            plt.ylabel('$\Delta$PROseq ({} over {})'.format( cellType_labels[treatment_y],cellType_labels[control_y]))
            # plt.title('$\Delta${}'.format(factor))
            x_sort = np.sort(x)
            plt.plot(x_sort,x_sort*slope+intercept,c = 'k',ls='--',lw=.9)
            plt.text(.97,.97,'$r={:.2f}$'.format(r_value),fontsize=12,transform=plt.axes().transAxes,ha='right',va='top')
            # plt.text(.50,1.03,'$P={:.1e}$ '.format(p_value),fontsize=13,transform=plt.axes().transAxes)
            plt.savefig(outdir+os.sep+'{}.pdf'.format(output_prename),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
    
    
# proseq_RPKM_df.to_csv('')
            
    
