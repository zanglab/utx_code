import sys,argparse
import os,glob
import numpy as np
import pandas as pd


def return_diff_genes_deseq2(csv_file,adjp,logfc):

    deseq_out = pd.read_csv(csv_file)
    deseq_out_upgenes = deseq_out.loc[(deseq_out['log2FoldChange']>logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_upgenes = set(deseq_out_upgenes['GeneID'])

    deseq_out_dngenes = deseq_out.loc[(deseq_out['log2FoldChange']<-1*logfc) & (deseq_out['padj']<adjp)]#.tolist();
    deseq_out_dngenes = set(deseq_out_dngenes['GeneID'])  
    return deseq_out_upgenes,deseq_out_dngenes


def write_DEG(genes,outfile):
    with open(outfile,'w') as outf:
        outf.write('\n'.join(genes));
        if len(genes)>0:
            outf.write('\n') # shell use \n to count the lines in each file
                
def main(indir,outdir):
    
    os.makedirs(outdir,exist_ok=True)

    deseq2_files = glob.glob(indir+'/*.csv')
    for deseq2_file in deseq2_files:
        basename = os.path.basename(deseq2_file).split('.txt')[0]
        
        adjp,logfc = 0.01,1
        upgenes,dngenes = return_diff_genes_deseq2(deseq2_file,adjp,logfc)
        write_DEG(upgenes,outdir+os.sep+basename+'_adjp{}_logfc{}_upgenes.txt'.format(adjp,logfc))
        write_DEG(dngenes,outdir+os.sep+basename+'_adjp{}_logfc{}_dngenes.txt'.format(adjp,logfc))
            
        adjp,logfc = 0.05,0.58
        upgenes,dngenes = return_diff_genes_deseq2(deseq2_file,adjp,logfc)
        write_DEG(upgenes,outdir+os.sep+basename+'_adjp{}_logfc{}_upgenes.txt'.format(adjp,logfc))
        write_DEG(dngenes,outdir+os.sep+basename+'_adjp{}_logfc{}_dngenes.txt'.format(adjp,logfc))
            
        #print(upgenes,basename);exit(0)



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
