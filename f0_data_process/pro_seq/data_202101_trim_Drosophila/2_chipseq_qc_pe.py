import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
from get_reads_positions import reads_positions
from IOparser_BedBam import get_tag_regions
import subprocess
import time
import get_RPM_RPKM_count_on_regions
import find_overlap_keep_info_NOT_sep_strand_asimport

#import re,bisect
plus = re.compile('\+')
#minus = re.compile('\-')
name = re.compile('\#name')
treat = re.compile('\#treat')
control = re.compile('\#control')

def get_lines(infile):

    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines


def write_fq_process_qc(outname,splitfile1,fq_base,fqdir,outdir,species):
    # write the process of each of the fq file
    fqfile = fqdir+os.sep+splitfile1;print(fqfile)
    bamfile = outdir+os.sep+outname+'_'+fq_base+'.bam'
    #bedfile = outdir+os.sep+fq_base+'.bed'
    time1 = time.time()
    if fqfile.split('.')[-1]=='gz':
        total_reads = int(subprocess.Popen('zcat {}|wc -l'.format(fqfile),shell=True,stdout=subprocess.PIPE).stdout.read().strip())/4
    elif fqfile.split('.')[-1]=='fastq':
        total_reads = get_lines('{}'.format(fqfile))/4
    time2 = time.time();print('fq',time2-time1)
    # num of mapped reads in bam file
    mapped_reads = int(subprocess.Popen('samtools view {}|wc -l'.format(bamfile),shell=True,stdout=subprocess.PIPE).stdout.read().strip())
    # num of uniq mapped reads
    time3 = time.time();print('bam',time3-time2)
    regions1, regions2 = get_tag_regions(species,'bam',bamfile)    
    uniq_reads=0
    for chrom in regions1.keys():
        for inner in regions1[chrom].keys():
            uniq_reads+=len(regions1[chrom][inner])
    for chrom in regions2.keys():
        for inner in regions2[chrom].keys():
            uniq_reads+=len(regions2[chrom][inner])
    time4 = time.time();print('uniq',time4-time3)
    return int(total_reads),int(mapped_reads),int(uniq_reads)



def write_out_qc(outf,block,fqdir,outdir,species,blockname,split_flag1,split_flag2):
    
    block = [i for i in block.split('\n') if len(i)>0]
    controlline = len(block)
    if len(block)>0:
        #print(block)
        for i in np.arange(len(block)):
            if name.match(block[i]):
                nameline = i
            if treat.match(block[i]):
                treatline = i
            if control.match(block[i]):
                controlline = i
        outname = block[nameline+1]
        treatnames = block[treatline+1:controlline]
        controlnames = block[controlline+1:len(block)]
        outdir = outdir+os.sep+outname
        os.makedirs(outdir,exist_ok=True)
       
        match1 = re.compile(split_flag1)   
        match2 = re.compile(split_flag2)
            
        # if is the block for match slurm file
        if outname==blockname:
            # summary of treat files
            total = 0
            for fqname in treatnames:
                if match1.search(fqname) and not match2.search(fqname):
                    fq_base = match1.split(os.path.basename(fqname))[0] # basename for matched paired end datasets
                    splitfile1 = [i for i in treatnames if (match1.search(i) and re.compile(fq_base).search(i))][0]
                    splitfile2 = [i for i in treatnames if (match2.search(i) and re.compile(fq_base).search(i))][0]
                    total_reads,mapped_reads,uniq_reads = write_fq_process_qc(outname,splitfile1,fq_base,fqdir,outdir,species)
                    outf.write('\n{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'.format(fq_base,total_reads,mapped_reads,uniq_reads,mapped_reads/total_reads,uniq_reads/mapped_reads))
                    total+=uniq_reads
            outf.write('\t{}'.format(total))
            
            # num of peak file    
            peak_file = '{}/{}_peaks.narrowPeak'.format(outdir,outname)
            peaks = get_lines(peak_file)
            outf.write('\t{}'.format(peaks))            
            # FRiP score
            merged_bam_file = '{}/{}_treat.bam'.format(outdir,outname)#;print(merged_bam_file)
            counting,bam_total = get_RPM_RPKM_count_on_regions.read_count_on_mapfile(peak_file,merged_bam_file,species,'bam')
            rip = sum([counting[ii][2] for ii in counting.keys()])
            outf.write('\t{:.3f}'.format(rip/bam_total)) 
            # peaks and peaks overlap UDHS
            udhs_file='/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'.format(species)
            peaks_overlapping_UDHS,nonoverlapped = find_overlap_keep_info_NOT_sep_strand_asimport.return_overlapped_count(peak_file,udhs_file,species)
            outf.write('\t{:.3f}'.format(peaks_overlapping_UDHS/max(peaks,1)))   # in case peaks =0

            # summary of control files            
            if len(controlnames)>0:
                total=0
                for fqname in controlnames:
                    if match1.search(fqname) and not match2.search(fqname):
                        fq_base = match1.split(os.path.basename(fqname))[0] # basename for matched paired end datasets
                        splitfile1 = [i for i in controlnames if (match1.search(i) and re.compile(fq_base).search(i))][0]
                        splitfile2 = [i for i in controlnames if (match2.search(i) and re.compile(fq_base).search(i))][0]
                        total_reads,mapped_reads,uniq_reads = write_fq_process_qc(outname,splitfile1,fq_base,fqdir,outdir,species)
                        outf.write('\n{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'.format(fq_base,total_reads,mapped_reads,uniq_reads,mapped_reads/total_reads,uniq_reads/mapped_reads))
                        total+=uniq_reads
                outf.write('\t{}'.format(total))



def main(infile,qc_out,species,blockname,split_flag1,split_flag2):
    
    outf = open(qc_out,'w')
    outf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('fastq file','total reads','mapped reads','non-redundant reads','mapping rate','non-redundant rate','total non-redundant reads','MACS2 peaks','FRiP','%peaks_on_UDHS'))            

    with open(infile) as inf:
        lines = inf.readlines()
        filedir = lines[0].strip().split('\t')[1]
        rundir = lines[1].strip().split('\t')[1]
        names_block = ''
        for line in lines[2:]:                   
            if len(line)>0:               
                if plus.match(line):
                    write_out_qc(outf,names_block,filedir,rundir,species,blockname,split_flag1,split_flag2)
                    names_block = ''                    
                else:
                    names_block = names_block+line
        write_out_qc(outf,names_block,filedir,rundir,species,blockname,split_flag1,split_flag2) 
    
    outf.write('\n')
    outf.close()    
            
   	        




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of match name list', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of slurm files,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    parser.add_argument('-b','--blockname', action = 'store', type = str,dest = 'blockname', help = 'block name for each slurm', metavar = '<str>')
    parser.add_argument('-f1', '--flag1', action = 'store', type = str,dest = 'flag1', help = 'split flag1 of input fastq files', metavar = '<file>')
    parser.add_argument('-f2', '--flag2', action = 'store', type = str,dest = 'flag2', help = 'split flag2 of input fastq files', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<3:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.outdir,args.species,args.blockname,args.flag1,args.flag2)
