import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def scatter_plot(df,marked_index,figname,mark_col):
    
    pylist = df.index
    plt.figure(figsize=(2.6,2.6))
    for index in df.index:
        xp = list(df.index).index(index)
        values = -1*np.log10(df.loc[index,mark_col])#;print(values)
        plt.scatter(xp,values,color='k',s=6)
    
    max_p = -1*np.log10(df.iloc[0,-1])
    values_reset = max_p*1.05
    # top 3 and TFs in shared top 10
    sorted_marker_ids = np.append(df.index[:10],df.index[:20].intersection(marked_index)[:5]);print(sorted_marker_ids)
    sorted_marker_pos = sorted(set([list(df.index).index(marker) for marker in  sorted_marker_ids]))
    x_pos_m=0
    for marker_id in sorted_marker_pos:
        xp = marker_id #;print(xp,df)
        values = -1*np.log10(df.loc[df.index[marker_id],mark_col])#;print(values)
        plt.scatter(xp,values,color='r',s=15)
        # ==== mark the index label, not overlap with each other
        values_reset = min(values,values_reset-max_p*0.095)#;print(values,values_reset)
        # ==== mark the index and plot the arrow separately
        if df.index[marker_id] in marked_index:
            plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id],color='red')
        else:
            plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id])
        
        # ==== plt.arrow(x,y,dx,dy)
        plt.arrow(xp+120+x_pos_m,values_reset+max_p*0.03,-90-x_pos_m,values-values_reset-max_p*0.03,\
                  length_includes_head = True,head_width=max_p*0.02,head_length=35,fc='k',ec='k')
        x_pos_m = x_pos_m + 30
    
    title = os.path.basename(figname).split('_bart3d_results.pdf')[0]
    title = title.split('_')
    plt.title('{}\n{}'.format(' '.join(title[:4]),' '.join(title[4:][::-1])),fontsize=12)
    plt.xlabel('TR Rank')
    plt.ylabel('-log$_{{10}}$ $p$-value')
    plt.axes().set_xticks([0,len(df.index)])
    plt.axes().set_xticklabels([1,len(df.index)],rotation=0, ha='center',fontsize=13,color='k')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()




def main_plot(infile,outdir):

    
    marked_index = []
    basename = os.path.basename(infile).split('.txt')[0]
    df = pd.read_csv(infile,index_col=0,sep='\t')
    mark_col = df.columns[-1]
    df = df.sort_values(by=[mark_col],ascending=True)
    figname = outdir+os.sep+basename+'.pdf'
    scatter_plot(df,marked_index,figname,mark_col)

    
def main(indir,outdir):
    
    os.makedirs(outdir,exist_ok=True)

    infiles = glob.glob('{}/*bart*_results.txt'.format(indir))
    for infile in infiles[:]:
        main_plot(infile,outdir)   




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
#     parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='rank_dot_figs')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outdir)
