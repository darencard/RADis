#!/usr/bin/env bash

usage()
{
cat << EOF
extract_fragments.sh

Version 1.0 (5 November, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script performs batch extraction of restriction enzyme fragments from a directory
containing a series of BED files (with ".bed" extension) and (if desired) performs a 
size selection. Output is both a BED file for the fragements and a FASTA file with the 
fragment sequences. The user must specify an input directory ('-i'). A prefix for the output
files can be specified using the '-o' flag. The lower and upper fragment range for size
selection is set using the '-l' and '-u' flags respectively (optional).

User must also specify the genome used for the restriction digest, which must have a
corresponding '.genomes' file for proper BED arrangement (see BEDtools documentation).

zsh perform_all_extractions.sh -g <genome> -i <input_dir> [-l <###> -u <###> -o <output> -h]

OPTIONS:
        -h      usage information and help (this message)
        -g      file name for the genome sequence in FASTA format (accompanying '.genomes'
                file must also be present)
	-i	the input directory containing the .bed files that are being summarized [.]	
	-l	the minimum (lower) fragment length to keep during size selection [NA]
	-u	the maximum (upper) fragment length to keep during size selection [NA]
        -o      directory for output files (files are named automatically) [output]
EOF
}

GENOME=
INPUT="."
LOWER=1
UPPER=10000000
OUTPUT="output"

while getopts "hg:i:l:u:o:" OPTION
do
        case $OPTION in
                h)
                        usage
                        exit 1
                        ;;
		g)
			GENOME=$OPTARG
			;;
		i)
			INPUT=$OPTARG
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

if [[ -z $GENOME ]] || [[ -z $INPUT ]]
then
	usage
	exit 1
fi

if [ ! -d "$OUTPUT" ]; then
	mkdir $OUTPUT
fi

# determine the number of bed files in the input directory and set
total=$(ls -1 $INPUT/*_[0-9].bed | wc -l)
# set total - 1 for outer loop
end1=$((total-1))

# for each file from 1 to n-1
for i in $(seq 1 $end1)
do
	# initialize second loop at first loop + 1
	j=$((i+1))
	# for each file from 2 to n
	for k in $(seq $j $total)
	do
		# set file 1 of comparison
		file1=$(ls $INPUT/*_[0-9].bed(.[$i]))
		# set core file name
		core=$(echo $file1 | rev | cut -d "_" -f 3- | rev)
		# gather the restriction enzyme name
		handle1=$(echo $file1 | rev | cut -d "_" -f 2 | rev)
		# gather the restriction enzyme length
		size1=$(echo $file1 | rev | cut -d "_" -f 1 | rev | cut -d "." -f 1)

		# set file 2 of comparison
		file2=$(ls $INPUT/*_[0-9].bed(.[$k]))
		# gather the restriction enzyme name
		handle2=$(echo $file2 | rev |cut -d "_" -f 2 | rev)
		# gather the restriction enzyme length
		size2=$(echo $file2 | rev | cut -d "_" -f 1 | rev | cut -d "." -f 1)

		# if length of file 1 restriction enzyme is greater than file 2 restrction enzyme
		if (($size1 >= $size2)); then
			# run with file 1 as rare and file 2 as common cutters
			zsh ./extract_fragments.sh -g $GENOME \
			-r $file1 -c $file2 -o $OUTPUT/$core"_"$handle1"-"$handle2
		# else if length of file 2 RE is greater than file 1 RE
		else;
			# run with file 2 as rare and file 1 as common cutters
			zsh ./extract_fragments.sh -g $GENOME \
			-r $file2 -c $file1 -o $OUTPUT/$core"_"$handle2"-"$handle1
		fi		
	done
done
