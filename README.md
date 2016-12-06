# RADis
RADis is a pipeline for performing in-silico RADseq experiments. It offers several advantages over other, existing software created for this purpose:

1. It performs all digests based on an input file of restriction enzymes of interest, and does not need to be run manually when evaluating multiple combinations of enzymes.
2. It exports information on the fragment coordinates and the fragment sequences in standard genomic data formats (BED and FASTA, respectively) that can be manipulated by existing powerful tools (e.g. BEDtools, etc.).
3. It is portable across most architectures, is easy to run, and works relatively quickly.

## Installation
RADis is composed of Python and Shell scripts, and should therefore work across most systems. The Shell scripts rely on the Z shell (zsh), so that must be available on the system. All Python code should be compatible with either Python 2 or Python 3.

[Heng Li's version of GNU sort](https://github.com/lh3/foreign/tree/master/sort), with support for alphanumerical sorting of genome coordinate files, is distributed with this software. I've modified the Makefile so that the executable name does not conflict with the existing version of sort. You must compile this software for certain RADis software to work.

```shell
make
```

Users must also temporarily or permanently add the executables to their system `$PATH`.
```shell
# copy to an existing path directory
cp restriction_digest_insilico.py extract_fragments.sh perform_all_extractions.sh sort_lh3 /existing/path/directory
# or export directory storing software to the $PATH
export PATH=$PATH:/path/to/RADis
```

## Workflow
RADis contains two core scripts for performing *in-silico* RADseq experiments. 

1. The first is a Python script that performs the restriction digest based on a provided genome file (in FASTA format) and a tab-delimited text file of restriction enzymes.

The restriction enzymes file must be provided in the following format:
```shell
<enzyme_name><tab><recognition_sequence><tab><cut_position>
```

The sequence os oriented from 5' to 3' and the cut position represents the last position of the recognition sequence before the enzyme cuts in the 5' to 3' direction (1-indexed). Here is an example with a commonly-used restriction enzyme, SbfI:
```shell
# note that comment lines with a '#' prefix are ignored
# 5'... CCTGCA|GG ...3'
# 3'... GG|ACGTCC ...5'
sbfI	CCTGCAGG	6
```

RADis is distributed with a pretty exhaustive list of restriction enzymes, which can be found in the `restrction_enzymes.txt` file.

restriction_digest_insilico.py can then be run on the genome of interest, as follows:
```shell
restriction_digest_insilico.py --input <genome.fasta> --enzymes <enzymes.txt> --output <output_prefix>
```

For each restriction enzyme, the string `_<enzyme_name>_<enzyme_seq>_<enzyme_seq_length>.bed` is added to the designed `prefix` to create the output file name. The output is in standard BED format.

2. From here, users can run the Shell script `extract_fragments.sh` to extract the coordinates and fragment sequences for each double digest combination of restriction enzymes, as follows:
```shell
extract_fragments.sh -g <genome.fasta> -r <rare_results> -c <common_results> -l <lower_length> -u <upper_length> -o <output_prefix>
```

`rare_results` and `common_results` represent the output BED files from `restriction_digest_insilico.py` for the rare cutting enzyme and common cutting enzyme, respectively. `lower_length` and `upper_length` indicate the window between which fragments of certain sizes will be kept (i.e., performs size selection'.

The output includes the coordinates of each double digest fragment in the genome in BED format and the resulting fragment sequences in FASTA format, both gzipped to save space.

Finally, a third script called `perform_all_extractions.sh` is also provided, allowing users to automate an analysis of all restriction enzyme pair combinations based on a directory of BED digest files from `restriction_digest_insilico.py`.
```shell
perform_all_extractions.sh -g <genome.fasta> -i <input_directory> -l <lower_length> -u <upper_length> -o <output_directory>
```

This script uses the standardized file names from `restriction_digest_insilico.py` to determine the rare and common cutters in each enzyme combination and runs `extract_fragments.sh` on all pairs.

## Features to Implement
1. Need to provide ability to parallelize restriction digests across multiple cores (right now it is defaulted to 4, but need to make code more flexible to choose).
2. Currently only works with double digest, so need to extend to create coordinates and fragments for a single digest.
3. Merge all functionality into one program with a configuration file specifying genomes and enzymes.
