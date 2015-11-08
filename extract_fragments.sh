#!/usr/bin/env bash

usage()
{
cat << EOF
extract_fragments.sh

Version 1.0 (1 November, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or 
desirable functioning.

This script extracts fragments from two restriction enzyme digest BED files and
(if desired) performs a size selection. Output is both a BED file for the fragements
and a FASTA file with the fragment sequences. The user should specify the rare cutting
enzyme using '-r' and the common-cutting enzyme using '-c'. A prefix for the output 
files can be specified using the '-o' flag. The lower and upper fragment range for size
selection is set using the '-l' and '-u' flags respectively (optional).

User must also specify the genome used for the restriction digest, which must have a
corresponding '.genomes' file for proper BED arrangement (see BEDtools documentation).

bash extract_fragments.sh -g <genome> -r <rare.bed> -c <common.bed> [-l <###> -u <###>
-o <output> -h]

OPTIONS:
	-h	usage information and help (this message)
	-g	file name for the genome sequence in FASTA format (accompanying '.genomes' 
		file must also be present)
	-r	BED file of restriction cut sites for the rare cutter (longer recognition 
		site)
	-c	BED file of restriction cut sites for the common cutter (shorter longer 
		recognition site)
	-l	the minimum (lower) fragment length to keep during size selection [NA]
	-u	the maximum (upper) fragment length to keep during size selection [NA]
	-o	output file name [output]
EOF
}

GENOME=
RARE=
COMMON=
LOWER=1
UPPER=10000000
OUTPUT="output"

while getopts "hg:r:c:l:u:o:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		g)
			GENOME=$OPTARG
			;;
		r)
			RARE=$OPTARG
			;;
		c)
			COMMON=$OPTARG
			;;
		l)
			LOWER=$OPTARG
			;;
		u)
			UPPER=$OPTARG
			;;
		o)
			OUTPUT=$OPTARG
			;;
		?)
			usage
			exit
			;;
	esac
done

if [[ -z $GENOME ]] || [[ -z $RARE ]] || [[ -z $COMMON ]]
then
	usage
	exit 1
fi

RENZ=$(echo $RARE | rev | cut -d "_" -f 3 | rev)
CENZ=$(echo $COMMON | rev | cut -d "_" -f 3 | rev)

# Find closest common cut site to each rare cut site using the respective BED files and format output into BED output
bedtools closest -g $GENOME.genomes -s -io -D a -t first -iu -fd -a $RARE -b $COMMON $RARE | \
grep -v "\-1" | awk '{if ($14 >= LOWER && $14 <= UPPER) print $0;}' LOWER="$LOWER" UPPER="$UPPER" - | \
awk '{if ($6 == "+") print $1 "\t" $2 "\t" $9 "\t" $4 "-" $11 "\t" $14 "\t" $6; \
else print $1 "\t" $9 "\t" $2 "\t" $4 "-" $11 "\t" $14 "\t" $6}' - | \
grep -v "$RENZ-$RENZ" - \
 > $OUTPUT.bed

# Use BED output from above and extract the fragment sequences in fasta format
bedtools getfasta -s -fi $GENOME -bed $OUTPUT.bed -fo $OUTPUT.fasta

# gzip both files to help save space
gzip $OUTPUT.bed
gzip $OUTPUT.fasta
