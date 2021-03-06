'''
Created on 04/15/2018
@authors: Chongzhi Zang, Zhenjia Wang<zw5j@virginia.edu>

This file is used to count reads from input tagsfile of bed/bam format 
for each start-end-region in input file
'''

import os,re,argparse,sys
import bisect
from IOparser_BedBam import get_tag_regions
#from BART.OptValidator import

plus = re.compile('\+')
minus = re.compile('\-')
      

def is_list_sorted(mylist):
    '''
    Check if list is sorted
    '''        
    for i in range(len(mylist)-1):
        if mylist[i] > mylist[i+1]:
            return 0
    return 1   
        	

def get_read_positions(positions,regions,val,fragment_size):
    '''
    Return all the shifted positions of reads 
    '''
    # default fragment_size = 150
    if fragment_size >= 0:
        shift = int(round(fragment_size/2))	
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                positions[chrom].extend([outer+shift*val for outer in regions[chrom][inner]])
    else:
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                for outer in regions[chrom][inner]:
                    shift = int(abs(outer-inner)/2)	
                    positions[chrom].append(outer+shift*val)
    return positions   
	        	    
		
def get_count_on_mapID(start,end,positions,mid,expand):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    if start < end and mid:
        middle = int((start+end)*0.5)
        s = bisect.bisect_left(positions,middle-expand)
        e = bisect.bisect_right(positions,middle+expand)
        return e-s, expand*2
        
    elif start < end:
        s = bisect.bisect_left(positions,start-expand)
        e = bisect.bisect_right(positions,end+expand)
        return e-s, end-start+expand*2
            
    else:
        return 0,1


def read_count_on_mapfile(infile,readsfile,species,format,fragmentsize=147,mid=False,expand=0):
    '''
    Count the num of (unique) reads in user-input bed/bam file
    on each of the UDHS -- bart profile
    '''
    #specify the species
    # get the start-end regions of each read in each chrom (as key) and separate by strand

    regions1, regions2 = get_tag_regions(species,format,readsfile)

    # get the counting position of each read(tag), chrom as key
    positions = get_read_positions({},regions1,1,fragmentsize)
    positions = get_read_positions(positions,regions2,-1,fragmentsize)
    total=0
    for chrom in positions:
        positions[chrom].sort() 
        total+=len(positions[chrom])
    # check if positions is NULL
    #if total ==0:
        #sys.stderr.write('Can not read the input bed/bam file!\n')
    # count reads on each DHS

    
    mapfile = open(infile,'r')
    line = mapfile.readline()
    counting = {}
    line_id=1
    while line:
        line = line.strip().split()
        #dhs = BED(line[0],line[1],line[2],line[3],line[4],line[5])
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        map_id = line[3]
#         try:
#             map_id = line[3]
#         except:
#             map_id = line_id
        if chrom in positions:
            #assert DHS_id not in counting
            nums,region_width = get_count_on_mapID(start,end,positions[chrom],mid,expand)
            rpkm=round(nums*1000000000/(total*region_width),3)
            rpm=round(nums*1000000/(total),3)
            counting[map_id]= [rpm,rpkm,nums]
        else:
            counting[map_id]=[0,0,0]
        line = mapfile.readline()
        line_id = line_id+1
    mapfile.close()
    return counting,total

def main(infile,readsfile,outfile,species,format,fragmentsize,mid,expand):

    counting,total = read_count_on_mapfile(infile,readsfile,species,format,fragmentsize,mid,expand)
    with open(outfile,'w') as outf:
        outf.write('ID\tRPM\tRPKM\tReadCount\n')
        for i in counting:
            outf.write('{}\t{}\t{}\t{}\n'.format(i,counting[i][0],counting[i][1],counting[i][2]))
    
    
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    parser.add_argument('-t', '--tagsfile', action = 'store', type = str,dest = 'tagsfile', help = 'ChIP-seq reads in bed/bam format', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile to write the the RPKM for each region', metavar = '<file>',required=True)
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
   
    parser.add_argument('-f', '--format', action = 'store', type = str,dest = 'format', help = 'format of file, bed or bam', metavar = '<str>',required=True)
    parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=False)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'expand of regions. Default: 0', metavar = '<int>',default=0)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.tagsfile,args.outfile,args.species,args.format,args.fragmentsize,args.mid,args.expand)
