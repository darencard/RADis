#!/usr/bin/env python

#print __name__

from __future__ import division
import optparse
import sys
from multiprocessing import Process
from Bio import SeqIO
from Bio import SeqUtils

usage_line = """
restriction_digest_insilico.py

Version 1.0 (1 November, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of \
proper or desirable functioning.

This script performs an in-silico restriction enzyme digest on an input sequence \
(i.e., a genome) using restriction enzymes supplied by the user in the form of a \
text file. The input sequence should be in FASTA format. The enzyme list is a \
tab-delimited file with 3 fields: (1) the enzyme name, (2) the enzyme sequence \
(5' to 3') and (3) the location in the enzyme sequence where the cut occurs \
(starting at 5' end, 0 indexed). Here is an example for a few restriction enzymes:
ecoRI	GAATTC	1
mspI	CCGG	1
sbfI	CCTGCAGG	6

Note: This script can handle ambiguity codes that are associated with some \
restriction endonucleases. It cannot handle lower-case letters, and will \
therefore miss restriction sites in soft-masked genomes.

The user can specify an output prefix that will be used for each output file. \
The output is a BED formatted file for each restriction enzyme. For each cut \
the BED will report both the sense and antisense cut sites.

python restriction_digest_insilico.py --input <sequence.fasta> --output <prefix> \
--enzymes <enzyme_list.txt>
"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", 
	action = "store", 
	type = "string", 
	dest = "input", 
	help = "input genome/sequence (in fasta format)")
parser.add_option("--output", 
	action = "store", 
	type = "string", 
	dest = "output", 
	help = "output file name prefix for gzipped tab-delimited coordinate results",
	default = "output")
parser.add_option("--enzymes", 
	action = "store", 
	type = "string", 
	dest = "enzymes", 
	help = "an input files with the restriction enzymes to use")

options, args = parser.parse_args()


##############################################################
### 	Parse restriction enzyme info from enzymes file    ###
##############################################################
		
def parse_RE_list(list, input):
	# for line in restriction enzyme input file
	for line in open(input, "r"):
		# only take lines that don't start with '#'
		if not line.strip().startswith("#"):
			# split by tab and append each list to overall restriction enzyme list
			bar = line.rstrip().split("\t")
			# 3 resulting fields: enzyme name, enzyme sequence, 5'-3' cut location (0 index)
			list.append(bar)
	# return filled list
	return list



##############################################################
### 	Perform digest using given sequence and enzyme     ###
##############################################################

def digest(enzyme, sequence, outfile, count):
	# search input sequence using enzyme sequence and return results to 'matches'
	matches = SeqUtils.nt_search(str(sequence.seq).upper(), enzyme[1])

	# for each of the items in results 'matches' list from 2nd item on (first item is match string)
	for match in matches[1:]:
		# create line for match on query stand
		line1 = sequence.id+"\t"+str(int(match)+int(enzyme[2]))+"\t"+str(int(match)+int(enzyme[2]))+"\t"+enzyme[0]+"\tcut-"+str(count)+"\t+\n"
		# look for reverse complement
		line2 = sequence.id+"\t"+str(int(match)+int(len(enzyme[1])-int(enzyme[2])))+"\t"+str(int(match)+int(len(enzyme[1])-int(enzyme[2])))+"\t"+enzyme[0]+"\tcut-"+str(count)+"\t-\n"

		# if cut site is past halfway point in enzyme, we should output antisense cut first to keep output BED sorted
		if len(enzyme[1])/2 < int(enzyme[2]):
			outfile.write(line2+line1)
		# if cut site is not past halfway point in enzyme, we can output in logical order
		else:
			# write both lines to ouput
			outfile.write(line1+line2)
		
		count += 1
	return count


##############################################################
###         Split Enzyme List for Parallelization          ###
##############################################################

def split_list(alist, wanted_parts=1):
	length = len(alist)
	return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
		for i in range(wanted_parts) ]



##############################################################
### 			Main Program		  	   ###
##############################################################

def main():
	# initialize empty restriction enzyme list
	RE_list = []
	
	# parse restriction enzyme list using function
	enzyme_list = parse_RE_list(RE_list, options.enzymes)

	batch_list = split_list(enzyme_list, wanted_parts=4)

	sys.stderr.write("Using target sequences from "+options.input+".\n")
		
	# for each enzyme in enzyme list and for each sequence, run digest function
	def batch1():
		for enzyme in batch_list[0]:
			sys.stderr.write("Running in silico digest using "+enzyme[0]+" ("+enzyme[1]+").\n")
			# file output format = <output_prefix>_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed
			output = open(options.output+"_"+enzyme[0]+"_"+enzyme[1]+"_"+str(len(enzyme[1]))+".bed", "w")
			count1 = 0
			for sequence in SeqIO.parse(open(options.input), "fasta"):
				new_count = digest(enzyme, sequence, output, count1)
			count1 += new_count
			output.close()

	def batch2():
		for enzyme in batch_list[1]:
			sys.stderr.write("Running in silico digest using "+enzyme[0]+" ("+enzyme[1]+").\n")
			# file output format = <output_prefix>_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed
			output = open(options.output+"_"+enzyme[0]+"_"+enzyme[1]+"_"+str(len(enzyme[1]))+".bed", "w")
			count2 = 0
			for sequence in SeqIO.parse(open(options.input), "fasta"):
				new_count = digest(enzyme, sequence, output, count2)
			output.close()

	def batch3():
		for enzyme in batch_list[2]:
			sys.stderr.write("Running in silico digest using "+enzyme[0]+" ("+enzyme[1]+").\n")
                        # file output format = <output_prefix>_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed
			output = open(options.output+"_"+enzyme[0]+"_"+enzyme[1]+"_"+str(len(enzyme[1]))+".bed", "w")
			count3 = 0
			for sequence in SeqIO.parse(open(options.input), "fasta"):
				new_count = digest(enzyme, sequence, output, count3)
			output.close()

	def batch4():
		for enzyme in batch_list[3]:
			sys.stderr.write("Running in silico digest using "+enzyme[0]+" ("+enzyme[1]+").\n")
			# file output format = <output_prefix>_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed
			output = open(options.output+"_"+enzyme[0]+"_"+enzyme[1]+"_"+str(len(enzyme[1]))+".bed", "w")
			count4 = 0
			for sequence in SeqIO.parse(open(options.input), "fasta"):
				new_count = digest(enzyme, sequence, output, count4)
			output.close()

	p1 = Process(target=batch1)
	p2 = Process(target=batch2)
	p3 = Process(target=batch3)
	p4 = Process(target=batch4)
	p1.start()
	p2.start()
	p3.start()
	p4.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()

### Run Main Program ###
main()
