import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *







outdir = 'f1_extract_data'
os.makedirs(outdir,exist_ok=True)

project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
project_dir='/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'


proseq_names = ['PCDH','DEL3','FL',]
norm_patterns=['RPKM','rawCount']
prenames = ['H3K4me1_overlapping_H3K27ac_overlapping_UTX','H3K27ac_H3K4me1_with_UTX']


for proseq_name in proseq_names[:]:
    for rep_id in ['1PRO','2PRO'][:]:
        for norm_pattern in norm_patterns[:]:
            # == get the PROseq binding pattern from two replicates
            rep_name='{}{}'.format(proseq_name,rep_id)
            proseq_binding_file = '{}/f6_proseq/data_0x2_MAPQ10/f4_pattern_at_PROseq_union_peak/{}_PROseq_peak_es2kb_bin200_{}.csv'.format(project_dir,rep_name,norm_pattern)
            proseq_binding_df = pd.read_csv(proseq_binding_file,sep='\t',index_col=0)
            for prename in prenames[:]:
                proseq_id_file = 'data/overlapped/PROseq_on_{}.bed'.format(prename)
                proseq_id_df = pd.read_csv(proseq_id_file,sep='\t',index_col=3,header=None)            
                # selected binding pattern
                df = proseq_binding_df.loc[proseq_id_df.index]
                df.to_csv('{}/{}{}_PROseq_peak_es2kb_bin200_{}_by_{}.csv'.format(outdir,proseq_name,rep_id,norm_pattern,prename),sep='\t')

