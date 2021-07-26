import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})


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


def compr_plot(x,y,lx,ly,thre_x,thre_y,elements,label_colors,figname,order=None,re_labels=None,mark_text=False):

    plt.figure(figsize=(4,4))
    label_mark=[]
    for ploc in np.arange(x.shape[0]):
        plabel = label_colors[elements[ploc]][1] if label_colors[elements[ploc]][1] not in label_mark else ''
        label_mark.append(label_colors[elements[ploc]][1])
        plt.scatter(x[ploc],y[ploc],c=label_colors[elements[ploc]][0],label=plabel,s=10)
        if mark_text:
            plt.axes().text(x[ploc],y[ploc],elements[ploc],fontsize=6)
    # ==== re-order the legends, if needed
    handles,labels = plt.axes().get_legend_handles_labels()
    if order:
        handles,labels = re_legend_labels(handles,labels,order,re_labels)
    plt.legend(handles, labels,bbox_to_anchor=[1,1],fontsize=16,frameon=False,\
               borderaxespad=0.1,labelspacing=.1,handletextpad=0.1,markerscale=2,ncol=1)
    if thre_x:
        plt.axes().axvline(thre_x,color='k',lw=1,ls='--')
    if thre_y:
        plt.axes().axhline(thre_y,color='k',lw=1,ls='--')
    plt.xlabel(lx)
    plt.ylabel(ly)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()


def get_label_color(column,label_colors_matchness):
    for key in label_colors_matchness.keys():
        if re.search(key,column) :
            return [label_colors_matchness[key],key]


def main(infile,outdir):
    
    os.makedirs(outdir,exist_ok=True)
    df = pd.read_csv(infile,index_col=0)
    label_colors_matchness = {'Chernobyl1245':'navy',\
                              'C1477':'skyblue',\
                              'C1575':'orange',\
                              'C1819':'red',\
                              'C1926':'purple'}
    # ==== initiate label color
    label_colors = {}
    for ele in df.index:
        label_colors[ele] = get_label_color(ele,label_colors_matchness)   

    order = [4,0,1,2,3]
    re_labels= {'Chernobyl1245':'diagnosis-1245',\
                'C1477':'relapse-tp1-1477',\
                'C1575':'relapse-tp2-1575',\
                'C1819':'relapse-tp3-1819',\
                'C1926':'relapse-tp4-1926'}


#     figname = '{}/qc_compr_uniqReads_vs_mappingRate.pdf'.format(outdir)
#     x,y = df['%mapping'],df['uniq reads']/1000000,
#     lx,ly = 'mapping rate(%)','non-reduncant reads (M)',
#     thre_x,thre_y = .3,.1,
#     compr_plot(x,y,lx,ly,thre_x,thre_y,df.index,label_colors,figname,order,re_labels)

    figname = '{}/qc_compr_mappedReads_vs_mappingRate.pdf'.format(outdir)
    x,y = df['%mapping'],df['mapped reads']/1000000,
    lx,ly = 'mapping rate(%)','mapped reads (M)',
    thre_x,thre_y = .3,.1,
    compr_plot(x,y,lx,ly,thre_x,thre_y,df.index,label_colors,figname,order,re_labels)


    figname = '{}/qc_compr_uniqRate_vs_mappingRate.pdf'.format(outdir)
    x,y = df['%mapping'],df['%unique'],
    lx,ly = 'mapping rate(%)','non-redundant rates (%)',
    thre_x,thre_y = .3,None,
    compr_plot(x,y,lx,ly,thre_x,thre_y,df.index,label_colors,figname,order,re_labels)






if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
#     parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of salmon results', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<5:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.outdir)

