import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
#sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree

def hi_plot(df,k,labels,label_colors,figname,metric,name_id):

    plt.figure(figsize=(2,10))
    Z = linkage(df,method = 'average',metric = metric)
    color_cut=Z[-(k - 1),2]
    dn = dendrogram(Z,labels = labels,color_threshold=color_cut,orientation='left',leaf_font_size=10)
    
    if label_colors:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            #print(lbl)
            lbl.set_color(label_colors[lbl.get_text()])
    
    plt.ylabel(name_id)
    plt.xlabel('{}'.format(metric))
    plt.savefig('{}_{}.pdf'.format(figname,metric),bbox_inches = 'tight',pad_inches = .1)
    plt.close()


def sns_clustering(df,figname):
    plt.figure(figsize = (5,40))
    sns.clustermap(df,metric = 'correlation')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = 0.1)
             

def get_label_color(column):
    if re.search('K27ac',column):
        return 'red'
    elif re.search('K4me3',column):
        return 'darkred'
    elif re.search('K4me1',column):
        return 'darkorange'
    elif re.search('K27me3',column):
        return 'k'
    elif re.search('MLL4',column):
        return 'blue'
    else:
        return 'deepskyblue'



# ==== main
    
outdir = 'f2_clustering_fig_by_antibody'
os.makedirs(outdir,exist_ok=True)
    
data_dir='RPKM_csv'
name_split='_RPKM'    
for name_id in ['promoter','genome_5kb']:
    csv_files = glob.glob('{}/*{}.csv'.format(data_dir,name_id))#;print(csv_files)
    csv_total = pd.DataFrame()
    for csv_file in csv_files:
        csv_df = pd.read_csv(csv_file,index_col=0,header=None,sep='\t')
        if len(csv_df.index) !=0:
            basename = os.path.basename(csv_file).split(name_split)[0]
            csv_df.columns = [basename]
            csv_total = pd.concat([csv_total,csv_df],axis=1)
    
    #df_sum = csv_total.sum()       
    #print(df_sum[df_sum==0])
    #sns_clustering(csv_total,'{}/{}.pdf'.format(outdir,'sns_heatmap'))
    print(csv_total.shape)
    
    label_colors = {}
    for column in csv_total.columns:
        label_colors[column] = get_label_color(column)
    #print(csv_total.columns)
    hi_plot(np.transpose(csv_total),1,csv_total.columns,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc',),'euclidean',name_id)
    hi_plot(np.transpose(csv_total),1,csv_total.columns,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc'),'correlation',name_id)
    hi_plot(np.transpose(csv_total),1,csv_total.columns,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc'),'cosine',name_id)



