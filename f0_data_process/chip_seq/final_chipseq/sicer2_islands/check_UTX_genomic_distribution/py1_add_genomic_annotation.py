import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")



exon_df = pd.read_csv('UTX_exons_overlapped.bed',sep='\t',index_col=3,header=None).loc[:,:3]
intron_df = pd.read_csv('UTX_introns_overlapped.bed',sep='\t',index_col=3,header=None).loc[:,:3]
promoter_df = pd.read_csv('UTX_promoter_overlapped.bed',sep='\t',index_col=3,header=None).loc[:,:3]

annotation_df = pd.read_csv('WT_UTX_with_Vector_control_peaks.narrowPeak',sep='\t',index_col=3,header=None).loc[:,:3]
annotation_df.columns = ['chr','start','end']
annotation_df['id']=annotation_df.index


annotation_df.loc[intron_df.index,'annotation']='Intron'
annotation_df.loc[exon_df.index,'annotation']='Exon'
annotation_df.loc[promoter_df.index,'annotation']='Promoter'
annotation_df['annotation'] = annotation_df['annotation'].fillna('Intergenic')

# #print(annotation_df);exit()
annotation_df.to_csv('UTX_binding_with_gene_annotation.csv',index=False)

a = annotation_df[annotation_df['annotation']=='Promoter'].shape[0]
b = annotation_df[annotation_df['annotation']=='Exon'].shape[0]
c = annotation_df[annotation_df['annotation']=='Intron'].shape[0]
d = annotation_df[annotation_df['annotation']=='Intergenic'].shape[0]

print('Promoter',a)
print('Exon',b)
print('Intron',c)
print('Intergenic',d)

## plot pie chart
fig = plt.figure(figsize=(2.6,2.6))
labels = 'Promoter', 'Exon', 'Intron', 'Intergenic'
sizes = [a,b,c,d]
explode = (0, 0.0, 0, 0)
plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90,counterclock=False)
plt.axes().axis('equal') 
plt.savefig('genomic_annotation.pdf',bbox_inches='tight',pad_inches=0.1,dpi=600)
plt.show()
plt.close()






