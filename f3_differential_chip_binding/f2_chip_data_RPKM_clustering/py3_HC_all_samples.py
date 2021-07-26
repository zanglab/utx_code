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

    plt.figure(figsize=(2,16))
    Z = linkage(df,method = 'average',metric = metric)
    color_cut=Z[-(k - 1),2];print(Z.shape,len(labels))
    # dn = dendrogram(Z,labels = labels,color_threshold=color_cut,orientation='left',leaf_font_size=10)
    dn = dendrogram(Z,labels = labels.to_list(),color_threshold=color_cut,orientation='left',leaf_font_size=13)
    
    if label_colors:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            #print(lbl)
            lbl.set_color(label_colors[lbl.get_text()])
    
    plt.ylabel(name_id)
    plt.xlabel('{}'.format(metric))
    plt.savefig('{}_{}.pdf'.format(figname,metric),bbox_inches = 'tight',pad_inches = .1)
#     plt.show()
    plt.close()


def sns_clustering(df,figname):
    plt.figure(figsize = (5,40))
    sns.clustermap(df,metric = 'correlation')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = 0.1)
             

def get_label_color(column):
    # if re.search('K27ac',column):
    #     return 'red'
    # elif re.search('K4me3',column):
    #     return 'darkred'
    # elif re.search('K4me1',column):
    #     return 'darkorange'
    # elif re.search('K27me3',column):
    #     return 'k'
    # elif re.search('MLL4',column):
    #     return 'blue'
    # else:
    #     return 'deepskyblue'

    if column.endswith('UTX') or column.endswith('HA') or column.endswith('UTX4'):
        return 'deepskyblue'
    elif column.endswith('MLL4') or column.endswith('MLL4SC') or column.endswith('MLL4GK'):
        return 'b'
    elif column.endswith('K27AC') or column.endswith('H3K27ac'):
        return 'red'
    elif re.search('K4me3',column):
        return 'darkred'
    elif column.endswith('K4me1') or column.endswith('K4me1_rep') or column.endswith('K4M1'):
        return 'darkorange'
    elif column.endswith('K27me3') or column.endswith('K27M3'):
        return 'k'
    elif column.endswith('P300'):
        return 'g'
    else:
        return 'purple'
    



# ==== main
    
outdir = 'f3_HC_clustering_fig'
os.makedirs(outdir,exist_ok=True)


for name_id in ['Promoter','Promoter_es10kb'][:]:
    rpkm_file= 'f1_RPKM_collection/{}_RPKM.csv'.format(name_id)
    with open(rpkm_file) as tpm_inf:
        df = pd.read_csv(tpm_inf,index_col=0)
#     df = df[(df>1).sum(axis=1)>(df.shape[1]*.05)];print('df dim\t',df.shape)
    
    label_colors = {}
    for column in df.columns:
        label_colors[column] = get_label_color(column)
    #print(df.columns)
    df = np.transpose(df)
    hi_plot(df,1,df.index,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc',),'euclidean',name_id)
    hi_plot(df,1,df.index,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc'),'correlation',name_id)
    hi_plot(df,1,df.index,label_colors, '{}/{}_{}'.format(outdir,name_id,'hc'),'cosine',name_id)



