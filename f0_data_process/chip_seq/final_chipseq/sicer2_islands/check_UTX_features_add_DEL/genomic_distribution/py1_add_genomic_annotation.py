import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")


for binding_type in ['WT','DEL_specific']:
    count_df = pd.DataFrame()
    exon_df = pd.read_csv('overlapped/{}_exons_overlapped.bed'.format(binding_type),sep='\t',index_col=3,header=None).loc[:,:3]
    intron_df = pd.read_csv('overlapped/{}_introns_overlapped.bed'.format(binding_type),sep='\t',index_col=3,header=None).loc[:,:3]
    promoter_df = pd.read_csv('overlapped/{}_promoter_overlapped.bed'.format(binding_type),sep='\t',index_col=3,header=None).loc[:,:3]
    
    annotation_df = pd.read_csv('../data/WT_peaks.bed',sep='\t',index_col=3,header=None).loc[:,:3]
    if binding_type == 'DEL_specific':
        annotation_df = pd.read_csv('../data/DEL_NOT_overlap_with_WT.bed',sep='\t',index_col=3,header=None).loc[:,:3]
    annotation_df.columns = ['chr','start','end']
    annotation_df['id']=annotation_df.index
    
    
    annotation_df.loc[intron_df.index,'annotation']='Intron'
    annotation_df.loc[exon_df.index,'annotation']='Exon'
    annotation_df.loc[promoter_df.index,'annotation']='Promoter'
    annotation_df['annotation'] = annotation_df['annotation'].fillna('Intergenic')
    
    # #print(annotation_df);exit()
    annotation_df.to_csv('{}_binding_with_gene_annotation.csv'.format(binding_type),index=False)
    
    a = annotation_df[annotation_df['annotation']=='Promoter'].shape[0]
    b = annotation_df[annotation_df['annotation']=='Exon'].shape[0]
    c = annotation_df[annotation_df['annotation']=='Intron'].shape[0]
    d = annotation_df[annotation_df['annotation']=='Intergenic'].shape[0]
    
    print('Promoter',a)
    print('Exon',b)
    print('Intron',c)
    print('Intergenic',d)
    count_df.loc['Promoter',binding_type] = a
    count_df.loc['Exon',binding_type] = b
    count_df.loc['Intron',binding_type] = c
    count_df.loc['Intergenic',binding_type] = d
    count_df.to_csv('{}_count.csv'.format(binding_type))
    
    ## plot pie chart
    fig = plt.figure(figsize=(2.6,2.6))
    labels = 'Promoter ({})'.format(a), 'Exon ({})'.format(b), 'Intron \n({})'.format(c), 'Intergenic \n({})'.format(d)
    sizes = [a,b,c,d]
    explode = (0, 0.0, 0, 0)
    plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=False, startangle=90,counterclock=False)
    plt.axes().axis('equal') 
    # plt.title('{}\n'.format(binding_type))
    plt.savefig('{}_genomic_annotation.pdf'.format(binding_type),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    
    
    
    
    
    
