import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from collections import Counter
    



def return_gene_loop_info(df_file,flag,hm,out_df):
    
    df = pd.read_csv(df_file)
    df = df.loc[df[['anchor1_promoter_Overlapped_genes','anchor2_promoter_Overlapped_genes']].dropna(how='all').index]
    # get cell_type info
    cell_type_info = df['cell_type'].values
    counter = Counter(cell_type_info)
    
    for key_ii in counter.keys():
        out_df.loc[key_ii,'{}_{}'.format(hm,flag)] = counter[key_ii]






# ==== main()

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang'
outdir = 'f2_deg_diff_loop_info'
os.makedirs(outdir,exist_ok=True)


for hm in ['H3K4me3','H3K27ac']:
    out_df = pd.DataFrame()  
    
    loop_file = '../f1_diff_looping/f1_diff_looping/{}.csv'.format(hm)
    flag =  'all_gene'
    return_gene_loop_info(loop_file,flag,hm,out_df)
    
    loop_file = 'f1_loop_anchor_with_deg/{}_dn.csv'.format(hm)
    flag =  'dn_gene'
    return_gene_loop_info(loop_file,flag,hm,out_df)
    
    loop_file = 'f1_loop_anchor_with_deg/{}_up.csv'.format(hm)
    flag =  'up_gene'
    return_gene_loop_info(loop_file,flag,hm,out_df)
    
    # total interactions
    out_df.loc['total'] = out_df.sum()
	# selected patterns
    del_index = [i for i in out_df.index if re.search('DEL',i) and not re.search('WT',i)]
    wt_index = [i for i in out_df.index if not re.search('DEL',i) and re.search('WT',i)]
    del_wt_index = [i for i in out_df.index if re.search('DEL',i) and re.search('WT',i)]

    out_df.loc['DEL-only'] = out_df.loc[del_index].sum()
    out_df.loc['WT-only'] = out_df.loc[wt_index].sum()
    out_df.loc['DEL-WT'] = out_df.loc[del_wt_index].sum()
    
    for flag in ['all_gene','dn_gene','up_gene']:
        out_df['% {}_{}'.format(hm,flag)] = out_df['{}_{}'.format(hm,flag)]/out_df.loc['total','{}_{}'.format(hm,flag)]
	
    # save the file
    out_df.to_csv(outdir+os.sep+'{}_deg_diff_loop_info.csv'.format(hm))

