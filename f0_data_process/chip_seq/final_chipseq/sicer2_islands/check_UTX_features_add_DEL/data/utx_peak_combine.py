import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
from scipy import stats



def read_peak_file(peak_file):
    df = pd.read_csv(peak_file,sep='\t',index_col=3,header=None)
    df.index = [i.split('/')[-1] for i in df.index]
    df.columns = ['chr','start','end','score','strand','signalValue','pValue','qValue','peak']
    df.insert(3,'name',df.index)
    return df

    

writer = pd.ExcelWriter('summary_UTX_peaks_WT_DEL.xlsx')

df = read_peak_file('WT_peaks.bed')
df.to_excel(writer,'WT peaks',index=False)

df = read_peak_file('WT_NOT_overlap_with_DEL.bed')
df.to_excel(writer,'WT specific peaks',index=False)

df = read_peak_file('WT_overlap_with_DEL.bed.uniq')
df.to_excel(writer,'WT overlap with DEL peaks',index=False)

df = read_peak_file('DEL_peaks.bed')
df.to_excel(writer,'DEL peaks',index=False)

df = read_peak_file('DEL_NOT_overlap_with_WT.bed')
df.to_excel(writer,'DEL specific peaks',index=False)

df = read_peak_file('DEL_overlap_with_WT.bed.uniq')
df.to_excel(writer,'DEL overlap with WT peaks',index=False)

# ++++ ChIP-seq binding at UTX sites
infile = '/Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/f0_data_integration_v3_RPKM/f3_combined_data_RPKM/ChIP-seq_binding_at_WT_UTX_peaks.xlsx'
df = pd.read_excel(infile,index_col=0)
df.to_excel(writer,'ChIP-seq binding at WT UTX')

writer.save()

