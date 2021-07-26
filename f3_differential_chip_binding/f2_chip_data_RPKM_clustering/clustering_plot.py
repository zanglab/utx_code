import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})
from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def re_legend_labels(handles,labels,order,re_labels):
    # === re-order the legend, if needed -- start
#     order = [4,0,1,2,3]
#     re_labels= {'Chernobyl1245':'diagnosis-1245',\
#                 'C1477':'relapse-tp1-1477',\
#                 'C1575':'relapse-tp2-1575',\
#                 'C1819':'relapse-tp3-1819',\
#                 'C1926':'relapse-tp4-1926'}
    handles,labels = [handles[idx] for idx in order],[re_labels[labels[idx]] for idx in order]
    return handles,labels
    
    

def pca_with_tsne(df,columns,label_colors,figname,order=None,re_labels=None,mark_text=False):

    pca=PCA(n_components=df.shape[0])
    projected_pca = pca.fit_transform(df)#;print(projected);exit()
    ratios = pca.explained_variance_ratio_
    print('pca+tsne pca dim:\t',projected_pca.shape)
    # ==== select PCs capture >50% variance
    for ii in np.arange(df.shape[0]):
        if sum(ratios[:ii])>.5:
            X = projected_pca[:,:ii]
            break
    projected = TSNE(n_components=2).fit_transform(X)
    print('pca+tsne dim:\t',projected.shape)
    
    plt.figure(figsize=(4,4))
    label_mark=[]
    for ploc in np.arange(projected.shape[0]):
        plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        label_mark.append(label_colors[columns[ploc]][1])
        plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]][0],label=plabel,s=10)
        if mark_text:
            plt.axes().text(projected[ploc,0],projected[ploc,1],columns[ploc],fontsize=6)
    
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if order:
        handles,labels = re_legend_labels(handles,labels,order,re_labels)

    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
               borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('tSNE_1')
    plt.ylabel('tSNE_2')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()



def tsne_plot(df,columns,label_colors,figname,order=None,re_labels=None,mark_text=False):

    projected = TSNE(n_components=2).fit_transform(df)
    print('tsne dim:\t',projected.shape)
    
    plt.figure(figsize=(4,4))
    label_mark=[]
    for ploc in np.arange(projected.shape[0]):
        plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        label_mark.append(label_colors[columns[ploc]][1])
        plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]][0],label=plabel,s=10)
        if mark_text:
            plt.axes().text(projected[ploc,0],projected[ploc,1],columns[ploc],fontsize=6)
    
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if order:
        handles,labels = re_legend_labels(handles,labels,order,re_labels)

    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
               borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('tSNE_1')
    plt.ylabel('tSNE_2')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()
    
    


def pca_plot(df,columns,label_colors,figname,order=None,re_labels=None,mark_text=False):

    pca=PCA(n_components=2)
    projected = pca.fit_transform(df)#;print(projected);exit()
    ratios = pca.explained_variance_ratio_
    print('pca dim:\t',projected.shape)
    
    plt.figure(figsize=(4,4))
    label_mark=[]
    for ploc in np.arange(projected.shape[0]):
        plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        label_mark.append(label_colors[columns[ploc]][1])
        plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]][0],label=plabel,s=10)
        if mark_text:
            plt.axes().text(projected[ploc,0],projected[ploc,1],columns[ploc],fontsize=6)
    
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if order:
        handles,labels = re_legend_labels(handles,labels,order,re_labels)

    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
               borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('PC1 ({:.1f}%)'.format(ratios[0]*100))
    plt.ylabel('PC2 ({:.1f}%)'.format(ratios[1]*100))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()




def hi_plot(df,labels,label_colors,figname,k=1):

    plt.figure(figsize=(3,15))
    Z = linkage(df) # method = 'average/single', metric = 'euclidean/correlation/cosine'
    color_cut=Z[-(k - 1),2]
    dn = dendrogram(Z,labels = labels,color_threshold=color_cut,orientation='left',leaf_font_size=5)
    
    if label_colors:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            #print(lbl,label_colors[lbl.get_text()]);exit()
            lbl.set_color(label_colors[lbl.get_text()][0])
    
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()



def sns_clustering(df,label_colors,figname,var_top=40000):
    # select top genes and plot heatmap
    df = df.dropna()
#     df_var = df.var(axis=1).sort_values(ascending=False)
#     df_var_top = df.loc[(df_var.index)[:np.min(var_top,len(df.index))]]
    
    import seaborn as sns
    sns.set(font_scale=.8)
    sns.set_style("whitegrid", {'axes.grid' : False})
    
    plt.figure()
    g = sns.clustermap(df,yticklabels=False,cmap = plt.cm.PiYG_r)
    for tick_label in g.ax_heatmap.axes.get_xticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(label_colors[tick_text][0])

    plt.savefig(figname,bbox_inches = 'tight',pad_inches = 0.1)
    plt.close()
             

