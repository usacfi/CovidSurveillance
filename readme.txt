This program aligns multiple COVID DNA sequences, trims the genes of interest, translates to protein, and identifies mutations in epitope regions.

Download the entire DNASeqAnalysis folder

# Example commandline to run script:

>>> python3 main.py -i references/Sequences/Morocco/48_Morocco_gisaid_hcov-19_2020_07_21_03.fasta  -aln True

Required parameters:
 -i is the input fasta file with the reference genome and test genomes 

Optional parameters:
 -ref is the reference genome name (Wuhan strain) the default is the name of the first genome in the input fasta file.
 -aln TRUE when the input needs alignment and FALSE when the input is already aligned
 -eps is the input text file that lists the epitope regions
 -loc is the input text file that indicates what and where are the relevant reference genes
 -g fasta file with only the reference genome
 -o is the output fasta file

# The program should run with Python 3.6 or later

Input files: 
 - 48_Morocco_gisaid_hcov-19_2020_07_21_03.fasta
 - Reference_gene_locations.txt
 - NCBI Reference Sequence_NC_045512.2.fasta
 - epitopes.txt 

Output files:
 - 01_temp_aligned.fasta
 - 02_trimmed.fasta
 - 03_protein.fasta
 - 04_mutations.fasta
 - 05_aminoacid_replacements.csv
 - 06_mutation_profile.csv
 - 07_mutation_profile.pdf

Fasta files can be viewed using any alignment viewers, e.g. AliView

Possible run-time errors:

1. Running the program in Windows or Linux might yield an error when (-aln TRUE). This is because the current MUSCLE and MAFFT programs in the folder are only compatible for MAC computers. 
Solution: You can download the MUSCLE/MAFFT programs **compatible to your OS** free online. Or you can change the input parameter (-aln FALSE), but note here that the program assumes that your input file is already aligned.
	
2. An error might come up when you are running in Python 2.
Solution: install Python3 (recommended) or convert the codes into Python 2

3. An error might come up when you don't have the following python modules:
 - argparse
 - numpy
 - pandas
 - matplotlib
 - itertools
 - os
Solution: Install the modules using pip in the command line/ terminal: 
pip3 install [module_name]







Let me know if there are corrections.
-jonathan 03/22/2021