#!/usr/local/env python

#print __name__

from __future__ import division
import optparse
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

def digest(enzyme, sequence, outfile):
	# search input sequence using enzyme sequence and return results to 'matches'
	matches = SeqUtils.nt_search(str(sequence.seq), enzyme[1])
	
	# for each of the items in results 'matches' list from 2nd item on (first item is match string)
	for match in matches[1:]:
		# create line for match on query stand and also for reverse complement on alternate strand
		line1 = sequence.id+"\t"+str(int(match)+int(enzyme[2]))+"\t"+str(int(match)+int(enzyme[2]))+"\t"+enzyme[0]+"\t.\t+\n"
		line2 = sequence.id+"\t"+str(int(match)+int(len(enzyme[1])-int(enzyme[2])))+"\t"+str(int(match)+int(len(enzyme[1])-int(enzyme[2])))+"\t"+enzyme[0]+"\t.\t-\n"

		# if cut site is past halfway point in enzyme, we should output antisense cut first to keep output BED sorted
		if len(enzyme[1])/2 < int(enzyme[2]):
			outfile.write(line2+line1)
		# if cut site is not past halfway point in enzyme, we can output in logical order
		else:
			# write both lines to ouput
			outfile.write(line1+line2)



##############################################################
### 			Main Program		  	   ###
##############################################################

def main():
	# initialize empty restriction enzyme list
	RE_list = []
	
	# parse restriction enzyme list using function
	enzyme_list = parse_RE_list(RE_list, options.enzymes)
	
	# for each enzyme in enzyme list and for each sequence, run digest function
	for enzyme in enzyme_list:
		# file output format = <output_prefix>_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed
		output = open(options.output+"_"+enzyme[0]+"_"+enzyme[1]+"_"+str(len(enzyme[1]))+".bed", "w")
		for sequence in SeqIO.parse(open(options.input), "fasta"):
			digest(enzyme, sequence, output)
		output.close()

### Run Main Program ###
main()
