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

writer.save()

