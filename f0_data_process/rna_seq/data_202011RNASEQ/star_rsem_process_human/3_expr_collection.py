import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd


def main(indir,outdir):

    os.makedirs(outdir,exist_ok=True)
    mapping_dirs = glob.glob(indir+'/*')
    out_count = pd.DataFrame()
    out_tpm = pd.DataFrame()
    out_fpkm = pd.DataFrame()
    
    for mapping_dir in sorted(mapping_dirs):
        sample_id = os.path.basename(mapping_dir)
        print(mapping_dir,sample_id)
        # ==== mapping rate
        expr_file = mapping_dir+os.sep+'/{}_rsem.genes.results'.format(sample_id)
        with open(expr_file) as expr_f:
            expr_df = pd.read_csv(expr_f,sep='\t',index_col=0)
        out_count = pd.concat([out_count,expr_df['expected_count'].rename(sample_id)],axis=1)
        out_tpm = pd.concat([out_tpm,expr_df['TPM'].rename(sample_id)],axis=1)
        out_fpkm = pd.concat([out_fpkm,expr_df['FPKM'].rename(sample_id)],axis=1)

    out_count.to_csv(outdir+os.sep+'expected_count.csv')
    out_tpm.to_csv(outdir+os.sep+'tpm.csv')
    out_fpkm.to_csv(outdir+os.sep+'fpkm.csv')




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of salmon results', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<5:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outdir)
