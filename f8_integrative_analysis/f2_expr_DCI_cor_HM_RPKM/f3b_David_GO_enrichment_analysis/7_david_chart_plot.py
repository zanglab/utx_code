import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"


def chart_barplot_goterm(infile,outdir,fdr,maxnum):
    
    inf = pd.read_csv(infile,sep="\t",index_col = 1)
    basename = os.path.basename(infile).split('.')[0];
    inf = inf[inf["FDR"]<fdr]
    #go_filtered = [i for i in inf.index if re.search(r'GO|mmu',i) ]
    #inf = inf.loc[go_filtered]
    
    inf = inf.sort_values(by=['FDR'],ascending=True).iloc[:maxnum,:]
    #print(inf);exit(0)
    y = [-np.log10(i) for i in inf["FDR"]]
    if len(y)>0:
        print(basename)
        plt.figure()
        fig, ax = plt.subplots(figsize=(6,len(inf.index)*0.3))  
        plt.barh(np.arange(len(inf)),y,align='center', alpha=0.9,height=0.7,color='cornflowerblue')
        plt.yticks(np.arange(len(inf)),[i.split("~")[-1] for i in inf.index])
        #plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
        plt.xlim([0,max(y)*1.2])
        plt.ylim([-0.8,len(y)-0.2])
        plt.gca().invert_yaxis()
        #plt.title(infile,verticalalignment='bottom')
        plt.xlabel("- log 10 (FDR)",fontsize=16,verticalalignment='top')
        for i,v in enumerate(y):
            ax.text(v + 0.03*y[0],i +0.4,str('{:.1f}'.format(v)))
        plt.title(basename[6:-15],va='bottom')
        plt.savefig(outdir+os.sep+'{}_GoTerm.pdf'.format(basename),bbox_inches='tight',pad_inches=0.1)



def chart_barplot_pathway(infile,outdir,fdr,maxnum):
    
    inf = pd.read_csv(infile,sep="\t",index_col = 1)#;print(inf)#;exit(0)
    basename = os.path.basename(infile).split('.txt')[0];#print(basename)
    inf = inf[inf["FDR"]<fdr]
    go_filtered = [i for i in inf.index if re.search('mmu',i)]
    inf = inf.loc[go_filtered]
    
    inf = inf.sort_values(by=['FDR'],ascending=True).iloc[:maxnum,:]
    #print(inf);exit(0)
    y = [-np.log10(i) for i in inf["FDR"]]
    if len(y)>0:
        print(basename)
        plt.figure()
        fig, ax = plt.subplots(figsize=(6,len(inf.index)*0.3))  
        plt.barh(np.arange(len(inf)),y,align='center', alpha=0.9,height=0.7,color='cornflowerblue')
        plt.yticks(np.arange(len(inf)),[i.split(":")[-1] for i in inf.index])
        #plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
        plt.xlim([0,max(y)*1.2])
        plt.ylim([-0.8,len(y)-0.2])
        plt.gca().invert_yaxis()
        #plt.title(infile,verticalalignment='bottom')
        plt.xlabel("- log 10 (FDR)",fontsize=16,verticalalignment='top')
        for i,v in enumerate(y):
            ax.text(v + 0.03*y[0],i +0.4,str('{:.1f}'.format(v)))
        plt.title(basename)
        plt.savefig(outdir+os.sep+'{}_KeggPathway.pdf'.format(basename),bbox_inches='tight',pad_inches=0.1)



def main(indir,outdir):

    os.makedirs(outdir,exist_ok=True)
    
    infiles = glob.glob(indir+'/*txt')
    for infile in infiles:
        chart_barplot_goterm(infile,outdir,0.05,20)
        # chart_barplot_pathway(infile,outdir,0.05,20)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outdir)