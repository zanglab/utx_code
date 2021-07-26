import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})
# import clustering_plot
from collections import Counter
from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def re_legend_labels(handles,labels,re_order,re_labels):
    # === re-order the legend, if needed -- start
#     order = [4,0,1,2,3]
#     re_labels= {'Chernobyl1245':'diagnosis-1245',\
#                 'C1477':'relapse-tp1-1477',\
#                 'C1575':'relapse-tp2-1575',\
#                 'C1819':'relapse-tp3-1819',\
#                 'C1926':'relapse-tp4-1926'}
    handles,labels = [handles[idx] for idx in re_order],[re_labels[labels[idx]] for idx in re_order]
    return handles,labels


def pca_plot(df,columns,label_colors,figname,re_order=False,re_labels=False,mark_text=False):

    pca=PCA(n_components=2)
    projected = pca.fit_transform(df)#;print(projected);exit()
    ratios = pca.explained_variance_ratio_
    print('pca dim:\t',projected.shape)
    
    plt.figure(figsize=(4,4))
    label_mark=[]
    for ploc in np.arange(projected.shape[0]):
        #print(label_mark)
        plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        label_mark.append(label_colors[columns[ploc]][1])
        plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]][0],label=plabel,s=30)
        if mark_text:
            plt.axes().text(projected[ploc,0],projected[ploc,1],columns[ploc],fontsize=7)
    
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if re_order:
        handles,labels = re_legend_labels(handles,labels,re_order,re_labels)

    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=12,frameon=True,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,markerscale=1.2,ncol=1)
    plt.xlabel('PC1 ({:.1f}%)'.format(ratios[0]*100))
    plt.ylabel('PC2 ({:.1f}%)'.format(ratios[1]*100))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()


def pca_with_tsne(df,columns,label_colors,figname,re_order=False,re_labels=False,mark_text=False):

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
        plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]][0],label=plabel,s=16)
        if mark_text:
            plt.axes().text(projected[ploc,0],projected[ploc,1],columns[ploc],fontsize=7)
    
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if re_order:
        handles,labels = re_legend_labels(handles,labels,re_order,re_labels)

    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=12,frameon=True,\
               borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('tSNE_1')
    plt.ylabel('tSNE_2')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()



def log_quantile_normalization(df):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = np.log2(df+1)
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized
             
def col_row_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    df = df.subtract(df.mean(axis=1),axis=0)
    return df

def log_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    return df


def get_label_color(column,label_colors_matchness):
    for key in label_colors_matchness.keys():
        if re.search(key,column) :
            return [label_colors_matchness[key],key]


def main():
    
    outdir = 'f5_PCA_QN'
    os.makedirs(outdir,exist_ok=True)

    tpm_file= 'f3_expr/tpm.csv'
    with open(tpm_file) as tpm_inf:
        df = pd.read_csv(tpm_inf,index_col=0)
    print('df dim\t',df.shape)
    
#     outliers = outliers.split();print('#outliers:\t',len(outliers))
#     df = df[df.columns.difference(outliers)]
    
    label_colors_matchness = {'del_cIDR':'k',\
                              'del_TPR':'purple',\
                              'Vector':'b',\
                              'WT':'darkred',\
                              'UTX_FUSIDR':'g',\
                              'UTX_eIFIDR':'orange'}
    # ==== initiate label color
    label_colors = {}
    for column in df.columns:
        label_colors[column] = get_label_color(column,label_colors_matchness)   
    print(label_colors)
    # need to be expressed in at least 5% samples
#     df = df[(df>1).sum(axis=1)>(df.shape[1]*.05)];print('df dim\t',df.shape)
    df = df[df.min(axis=1)>1];print('df dim\t',df.shape)
#     df = np.log2(df+1)
    df = log_quantile_normalization(df)
    
    re_order = [2,3,5,4,1,0]
    re_labels= {'del_cIDR':'$\Delta$cIDR',\
                'del_TPR':'$\Delta$TPR',\
                'Vector':'Vector',\
                'WT':'WT',\
                'UTX_FUSIDR':'UTX-FUS$_{IDR}$',\
                'UTX_eIFIDR':'UTX-eIF$_{IDR}$'}
    #re_order=False 
             
#     figname = '{}/tpm_hierarchical.pdf'.format(outdir)
#     clustering_plot.hi_plot(np.transpose(df),df.columns,label_colors,figname)

    figname = '{}/tpm_PCA.pdf'.format(outdir)
    pca_plot(np.transpose(df),df.columns,label_colors,figname,re_order,re_labels,mark_text=False)

    figname = '{}/tpm_PCA_followed_tSNE.pdf'.format(outdir)
#     pca_with_tsne(np.transpose(df),df.columns,label_colors,figname,re_order,re_labels,mark_text=False)





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
   # parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
