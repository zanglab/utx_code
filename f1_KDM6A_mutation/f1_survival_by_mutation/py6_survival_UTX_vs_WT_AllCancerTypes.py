import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)

sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
# from stats_pvalues import irwin_hall_cdf
from lifelines.statistics import logrank_test
#from lifelines.plotting import add_at_risk_counts
from lifelines import KaplanMeierFitter


def survival_for_two(df,treat,legends,figname):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    
    t1 = df.loc[ix]['time']
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value;print('pvalue:\t{}'.format(pvalue))
    
    # survival curves
    plt.figure(figsize=(3.,3.))
    ax = plt.subplot(111)
    
    kmf_exp = KaplanMeierFitter()
    #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
    g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                    censor_styles={'ms': 9, 'marker': '+'},ci_show=False,c='k',ls='--')
    
    kmf_control = KaplanMeierFitter()
    #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
    g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                        censor_styles={'ms': 9, 'marker': '+'},ci_show=False,c='r',ls='-')
    
    handles, labels = ax.get_legend_handles_labels();print(labels)
    lg = ax.legend(handles[1::2], legends,loc='upper right',fontsize=14,frameon=False,borderaxespad=.1,handletextpad=.2,labelspacing=.2,handlelength=1.1)
    if pvalue<0.05:
         plt.axes().text(df['time'].max()*0.65,0.45,'p={:.1e}'.format(pvalue),fontsize=16,ha='center')
    plt.ylim([-0.02,1.05])
#     plt.xlim([0,max_val*1])
    plt.xlabel('Days',fontsize=22)
    plt.ylabel('Survival probability',fontsize=22)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
    plt.close()
    return results




def retuan_clinical_df(clinical_dir,project_code):
    # read clinical info and return formatted df
    df = pd.read_csv('../data/{}/donor.tsv'.format(clinical_dir),sep='\t',index_col=0)
#     df = df[df['project_code']==project_code]
    time_columns = ['donor_survival_time','donor_interval_of_last_followup']
    df['time'] = df[time_columns].max(axis=1)
    df['status'] = df['donor_vital_status']
    df.loc[df['status']=='deceased','status']=1
    df.loc[df['status']=='alive','status']=0
    df = df[['time','status']].dropna(how='any')
    df = df[df['time']>3];print(df.shape)
    return df
    



# == main ==

outdir = 'f6_figs_UTX_vs_WT_AllCancerTypes'
os.makedirs(outdir,exist_ok=True)

dir_all='icgc_clinical_donors_all_20200823'
dir_UTX='icgc_clinical_donors_with_mutation_KDM6A_20200823'
dir_UTX_cIDR='icgc_clinical_donors_with_mutation_KDM6A_cIDR_20200823'

project_code='BLCA-US'   

df_all = retuan_clinical_df(dir_all,project_code)
df_UTX = retuan_clinical_df(dir_UTX,project_code)
df_UTX_cIDR = retuan_clinical_df(dir_UTX_cIDR,project_code)

# == this is what to change
treat_ids = df_UTX.index
control_ids = df_all.index.difference(df_UTX.index) 
# control_ids = df_UTX.index.difference(df_UTX_cIDR.index) 
legends=['UTX not mutated','UTX mutated',]
print(len(treat_ids))
print(len(control_ids))

clinical_df = df_all.loc[treat_ids.union(control_ids)];print(clinical_df.shape)
clinical_df['group'] = 'control'
clinical_df.loc[treat_ids,'group']='treat'
clinical_df.loc[control_ids,'group']='control'

clinical_df.to_csv(outdir+os.sep+'survival_info.csv')
figname = outdir+os.sep+'survival_UTX_cIDR_vs_WT.png'
survival_for_two(clinical_df,'treat',legends,figname)


