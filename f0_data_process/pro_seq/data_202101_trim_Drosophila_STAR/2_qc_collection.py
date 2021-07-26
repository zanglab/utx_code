import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd


def main(indir,outfile):

    mapping_dirs = glob.glob(indir+'/*')
    outdf = pd.DataFrame(columns = ['total reads','mapped reads','%mapping','duplicates','%duplicates','uniq reads','%unique'])
    
    for mapping_dir in sorted(mapping_dirs):
        sample_id = os.path.basename(mapping_dir)
        print(mapping_dir,sample_id)
        # ==== mapping rate
        log_file = mapping_dir+os.sep+'/{}Log.final.out'.format(sample_id)
        if not os.path.isfile(log_file):
            print(log_file)
            continue
        with open(log_file) as log_f:
            lines = log_f.readlines()
        for line in lines:
            if re.search("Number of input reads",line):
                outdf.loc[sample_id,'total reads'] = float(line.strip().split()[-1])
            if re.search("Uniquely mapped reads number",line):
                outdf.loc[sample_id,'mapped reads'] = float(line.strip().split()[-1])

        # ==== non-duplicates
        log_file = mapping_dir+os.sep+'/{}Aligned.sortedByCoord.out.markdup.txt'.format(sample_id)
        with open(log_file) as log_f:
            lines = log_f.readlines()
        info = lines[7].strip().split('\t')
        outdf.loc[sample_id,'duplicates'] = float(info[5]) + float(info[6])
        outdf.loc[sample_id,'%duplicates'] = float(info[8])
        
    outdf['%mapping'] = outdf['mapped reads']/outdf['total reads']      
    outdf['uniq reads'] = outdf['mapped reads'] - outdf['duplicates']          
    outdf['%unique'] = outdf['uniq reads']/outdf['mapped reads']        
    outdf.to_csv(outfile)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of salmon results', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<5:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outfile)
