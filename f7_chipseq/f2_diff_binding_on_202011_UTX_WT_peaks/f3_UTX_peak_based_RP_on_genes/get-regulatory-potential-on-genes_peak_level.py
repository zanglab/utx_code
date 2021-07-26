#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from GenomeData import *
#import get_total_tag_counts
import SeparateByChrom
import associate_island_with_genes
import associate_tags_with_regions
import Utility

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq peak file")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-e", "--extension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>", default = 100000)
	parser.add_option("-f", "--decayfactor", action="store", type="int", dest="decayfactor", help="decay factor in RP calculation, default 4", metavar="<int>", default = 4)
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
		
	if Utility.fileExists(opt.bedfile):
		SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1');
	else:
		print opt.bedfile, " not found";
		sys.exit(1)
	
	#genomelength = 0.0
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	#total = 0;
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			(gene_name_list, region_start_list, region_end_list) = associate_island_with_genes.get_gene_feature_region_lists(coords[chrom], "Promoter", 0)
			
			tag_position_list = []
			tag_score_list = []
			read_file = chrom + ".bed1";
			f = open(read_file,'r')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_position_list.append(associate_tags_with_regions.tag_position(sline, -1))
					tag_score_list.append(atof(sline[8]))
			f.close();
			position_list = []
			score_list = []
			for (x,y) in sorted(zip(tag_position_list, tag_score_list)): 
				position_list.append(x)
				score_list.append(y)
			tag_count_list = associate_tags_with_regions.find_strength_distance_weighted_score_on_regions(position_list, score_list, region_start_list, opt.promoter_extension, opt.decayfactor)
			
			assert len(gene_name_list) == len(tag_count_list)
			for i in range(0, len(gene_name_list)):
				occupied_genes[gene_name_list[i]] = tag_count_list[i]
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()
	SeparateByChrom.cleanup(chroms, '.bed1');


if __name__ == "__main__":
	main(sys.argv)
