import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import AML_modules
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn2_plot(a,b,la,lb,figname,compr_type,deg_type,color1='purple',color2='skyblue'):
    # fisher exact test
    total=19826
    va = len(a.intersection(b))
    vb = len(a) - len(a.intersection(b))
    vc = len(b) - len(a.intersection(b))
    vd = total - len(a) - len(b) + len(a.intersection(b))
    s,p = stats.fisher_exact([[va,vb],[vc,vd]])
    print(compr_type,deg_type,len(a),len(b))
    print('odds ratio={:.2f}, pvalue={:.2e}'.format(s,p))
    
    # plot venn diagram
    if compr_type=='UTX-DEL3-72':
        color1='red'
    plt.figure(figsize=(4,4))
    out = venn2([a,b],set_labels=(la,lb),set_colors=(color1,color2), alpha=0.5)
    
    for text in out.set_labels:
        text.set_fontsize(16)
    for text in out.subset_labels:
        try:
            text.set_fontsize(16)
        except:
            pass
    if deg_type=='dn':
        deg_type = 'down'
    plt.title('{}-genes'.format(deg_type))
    if p<0.05:
        plt.text(x=.2,y=-.15,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches=0.1,transparent=True)
    plt.close()
    
    with open(figname+'.txt','w') as outf:
        outf.write('\n'.join(a.intersection(b))+'\n')




def main(infile):
    
    project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/'
    indir='{}/f0_data_process/rna_seq/data_202011RNASEQ/star_rsem_process_human/f6_deg/f2_deg'.format(project_dir)
    outdir = 'f2_human_rna_venn_DEG'
    os.makedirs(outdir,exist_ok=True)
    
    for deg_type in ['up','dn']:
        for compr_type in ['MT2_THP1']:
            deg_list_t = '{}/treated_{}_vs_ctrl_PCDH_THP1.deseq2.csv_adjp0.05_logfc0.32_{}genes.txt'.format(indir,compr_type,deg_type)
            deg_list_c = '{}/treated_UTXWT_THP1_vs_ctrl_PCDH_THP1.deseq2.csv_adjp0.05_logfc0.32_{}genes.txt'.format(indir,deg_type)
        
            deg_t = [i.strip() for i in open(deg_list_t).readlines()]
            deg_c = [i.strip() for i in open(deg_list_c).readlines()]
        
            figname = outdir+os.sep+'{}_{}_venn_DEG.png'.format(compr_type,deg_type)
            venn2_plot(set(deg_t),set(deg_c),compr_type,'WT',figname,compr_type,deg_type)
    






if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)
