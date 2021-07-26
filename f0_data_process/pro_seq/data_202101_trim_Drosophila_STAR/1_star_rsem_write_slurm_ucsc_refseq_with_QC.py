import os,sys,argparse
import fileinput,time
import glob,re
import re,bisect
import pandas as pd
import numpy as np
#def expand_region(summitlist):
#from reads_count import read_count_on_mapfile

#import re,bisect
plus = re.compile('\+')
#minus = re.compile('\-')
name = re.compile('\#name')
fastq_files = re.compile('\#files')


def write_slurm(block,slurmdir,filedir,salmondir,species):

    block = [i for i in block.split('\n') if len(i)>0]
    totallines = len(block)
    if len(block)>0:
        #print(block)
        for i in np.arange(len(block)):
            if name.match(block[i]):
                nameline = i
            if fastq_files.match(block[i]):
                fileline = i
        outname = block[nameline+1]
        filenames = block[fileline+1:totallines]
        filenames = [filedir+os.sep+'{}'.format(i) for i in filenames]
        #print(filenames)

        pair1_files = sorted([i for i in filenames if re.search('R1_001_val_1.fq.gz',i)])
        pair2_files = sorted([i for i in filenames if re.search('R2_001_val_2.fq.gz',i)])
        #print(pair2_files);exit(0)

        with open(slurmdir+os.sep+'{}.slurm'.format(outname),'w') as slurmout:
                slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
''')
                #######################
                # Need to change the index file and rsem reference
                #########################
              
                slurmout.write('#SBATCH -o {}/{}.out \n\n'.format(slurmdir,outname))
                index_file = '/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/{}'.format(species)
                rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/{}/rsem_reference".format(species)
                mapping_outdir = '{}/{}'.format(salmondir,outname)
                os.makedirs(mapping_outdir,exist_ok=True)
                slurmout.write('index_file="{}"\n'.format(index_file))
                slurmout.write('rsem_reference="{}"\n'.format(rsem_reference))
                slurmout.write('mapping_outdir="{}"\n\n'.format(mapping_outdir))
                # ==== STAR MAPPING
                slurmout.write('# ==== STAR MAPPING\n')
                slurmout.write('''STAR \\
--runMode alignReads \\
--genomeDir ${{index_file}} \\
--readFilesIn \\\n{} \\\n{} \\
--outFileNamePrefix ${{mapping_outdir}}/{} \\
--outFilterMultimapNmax 1 \\
--runThreadN 4 \\
--outSAMtype BAM SortedByCoordinate \\
--outSAMunmapped Within \\
--quantMode TranscriptomeSAM GeneCounts '''.format(',\\\n'.join(pair1_files),',\\\n'.join(pair2_files),outname))
                
                if pair1_files[0].endswith('.gz'):
                    slurmout.write('\\\n--readFilesCommand zcat\n\n')
                else:
                    slurmout.write('\n\n\n')
                
                # ==== mark duplicates
                slurmout.write('# ==== mark duplicates\n')
                slurmout.write('module load picard\n\n')
                slurmout.write('''java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \\
I=${{mapping_outdir}}/{}Aligned.sortedByCoord.out.bam \\
O=${{mapping_outdir}}/{}Aligned.sortedByCoord.out.markdup.bam \\
M=${{mapping_outdir}}/{}Aligned.sortedByCoord.out.markdup.txt \\
ASSUME_SORT_ORDER=coordinate \n\n\n'''.format(outname,outname,outname))


                # ==== RSEM quantification
                slurmout.write('# ==== RSEM quantification\n')
                slurmout.write('module load perl\n\n')
                slurmout.write('''rsem-calculate-expression \\
--num-threads 4 \\
--fragment-length-max 1000 \\
--estimate-rspd \\
--no-bam-output \\''')
                if len(pair2_files)>0:
                    slurmout.write('''
--paired-end \\
--bam ${{mapping_outdir}}/{}Aligned.toTranscriptome.out.bam \\
${{rsem_reference}} ${{mapping_outdir}}/{}_rsem \n\n'''.format(outname,outname))
                else:
                    slurmout.write('''
--bam ${{mapping_outdir}}/{}Aligned.toTranscriptome.out.bam \\
${{rsem_reference}} ${{mapping_outdir}}/{}_rsem \n\n'''.format(outname,outname))


            
def main(infile,slurmdir,species):

    os.makedirs(slurmdir,exist_ok=True)
    with open(infile) as inf:
        lines = inf.readlines()
        filedir = lines[0].strip().split('\t')[1]
        rundir = lines[1].strip().split('\t')[1]
        names_block = ''
        for line in lines[2:]:                   
            if len(line)>0:               
                if plus.match(line):
                    write_slurm(names_block,slurmdir,filedir,rundir,species)
                    names_block = ''                    
                else:
                    names_block = names_block+line
                    
        write_slurm(names_block,slurmdir,filedir,rundir,species)     


    
      
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of match name list', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of slurm files,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<3:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.outdir,args.species)
