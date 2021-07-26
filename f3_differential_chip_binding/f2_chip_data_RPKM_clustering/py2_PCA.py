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
# import clustering_plot
from collections import Counter
from sklearn.decomposition import PCA


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


# def get_label_color(column,label_colors_matchness):
#     for key in label_colors_matchness.keys():
#         if re.search(key,column) :
#             return [label_colors_matchness[key],key]



# 'data_1st_submit','data_20201127_1205_merged','data_20201209','data_20201220'

def get_label_marker(column):
    if re.search('data_1st_submit',column):
        return 'x'
    elif re.search('data_20201127_1205_merged',column):
        return 'o'
    elif re.search('data_20201209',column):
        return '^'
    elif re.search('data_20201220',column):
        return '*'
    
    
def get_label_color(column):
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
    
    
    
    
    
    
def pca_plot(df,columns,label_colors,label_markers,figname,order=None,re_labels=None,mark_text=True):

    pca=PCA(n_components=2)
    projected = pca.fit_transform(df)#;print(projected);exit()
    ratios = pca.explained_variance_ratio_
    print('pca dim:\t',projected.shape)
    
    plt.figure(figsize=(6,6))
    marker_legends,marker_list=[],[]
    color_legends=[]
    gs=[]
    for ploc in np.arange(projected.shape[0]):
        # plabel = label_colors[columns[ploc]][1] if label_colors[columns[ploc]][1] not in label_mark else ''
        # label_mark.append(label_colors[columns[ploc]][1])
        g = plt.scatter(projected[ploc,0],projected[ploc,1],c=label_colors[columns[ploc]],marker=label_markers[columns[ploc]], s=160,label=label_markers[columns[ploc]])
        # gs.append(g)
        # plt.scatter(projected[ploc,0],projected[ploc,1],s=10,marker='^')
        # if label_colors[columns[ploc]] not in marker_list:
        #     marker_legends.append(g)
        #     marker_list.append(label_colors[columns[ploc]])
        # if mark_text:
        text = columns[ploc].split('_')[-1]
        # plt.axes().text(projected[ploc,0],projected[ploc,1],text,fontsize=11)
    # g = plt.scatter(projected[:,0],projected[:,1],c=list(label_colors.values()),marker=list(label_markers.values()), s=110)
        
    # ==== re-order the legends, if needed
    # handles,labels = plt.axes().get_legend_handles_labels()
    # if order:
        # handles,labels = re_legend_labels(handles,labels,order,re_labels)
    # plt.legend(gs)
    # handles, labels = plt.axes().legend_elements(prop="sizes", alpha=0.6)
    # plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
                # borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    plt.xlabel('PC1 ({:.1f}%)'.format(ratios[0]*100))
    plt.ylabel('PC2 ({:.1f}%)'.format(ratios[1]*100))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
#     plt.show()
    plt.close()




# def main():


   
outdir = 'f2_PCA_figs'
os.makedirs(outdir,exist_ok=True)

for name_id in ['Promoter','Promoter_es10kb'][:]:
    rpkm_file= 'f1_RPKM_collection/{}_RPKM.csv'.format(name_id)
    with open(rpkm_file) as tpm_inf:
        df = pd.read_csv(tpm_inf,index_col=0)

    # ==== initiate label color
    label_colors = {}    
    label_markers = {}
    for column in df.columns:
        label_colors[column] = get_label_color(column)   
        label_markers[column] = get_label_marker(column)   

    # need to be expressed in at least 5% samples
    # df = df[(df>1).sum(axis=1)>(df.shape[1]*.05)];print('df dim\t',df.shape)
    # df = np.log2(df+1)
    df = log_quantile_normalization(df)
    # df = col_row_norm(df)
    # df = df.subtract(df.median(),axis=1)
    # df = df[[i for i in df.columns if label_colors[i]=='k']]
    # df = df[[i for i in df.columns if label_markers[i]!='x']]
    

#     figname = '{}/tpm_hierarchical.pdf'.format(outdir)
#     clustering_plot.hi_plot(np.transpose(df),df.columns,label_colors,figname)

    figname = '{}/{}_QN_PCA.pdf'.format(outdir,name_id)
    pca_plot(np.transpose(df),df.columns,label_colors,label_markers,figname)






# if __name__ == '__main__':
# 
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
#     #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
#     #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
#     #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
#    # parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
#     
# 
#     args = parser.parse_args()
#     if(len(sys.argv))<0:
#         parser.print_help()
#         sys.exit(1)
#   
#     main()
