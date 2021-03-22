# -*- coding:utf-8 -*-

import argparse
from classes.Sequence import Sequence
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from itertools import chain
import os

'''
Sample Run in Command Prompt/Terminal: 
python3 main.py -i references/Sequences/Morocco/48_Morocco_gisaid_hcov-19_2020_07_21_03.fasta -g input/reference_genes_locations.txt 
-o output/mutations.fasta -ref "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome" -aln True
'''


parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-i', '--input_alignments_fasta', dest='input', help='input alignment file', required=True)

parser.add_argument('-aln', '--needs_alignment', dest='needs_alignment', help='(Type: Bool) True if input still needs to be aligned,'
                    'False if input is already aligned, e.g. -aln True, default value is True', 
                    default=True, type=lambda x: (str(x).lower() in ['true', '1', 'yes']))

parser.add_argument('-loc', '--gene_loc', dest='gene_loc', help='input file with the location of genes', default='input/reference_genes_locations.txt' )

parser.add_argument('-o', '--output', dest='out', help='base name for output files', default='output/mutations.fasta')

parser.add_argument('-ref', '--reference_name', dest='ref_name',
                    help='name of reference sequence in alignment to compare with other sequences; otherwise uses first sequence in alignment as reference', 
                    default='NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome')

parser.add_argument('-g', '--genome', dest='ref_genome', help='fasta file with only the reference genome', 
                    default='references/Sequences/References\ Sequences/NCBI\ Reference\ Sequence_NC_045512.2.fasta')

parser.add_argument('-eps', '--epitope_regions', dest='eps', help='text file that lists epitopes to check', default='input/epitopes.txt')

args = parser.parse_args()



def parse_alignment(fasta_file, ref_name, gene_name):
  seq_dict = dict()

  with open(fasta_file) as data:
    name = ''
    seq = ''

    for line in data:
      line = line.rstrip() # remove trailing whitespaces

      if line.startswith('>'): # name
        if ref_name == None: # if user did not specify -ref input
          ref_name = line[1:] # assumes first name found is the reference
          name = gene_name
          seq = ''
          
        # if ref_name is defined by the user, look for it in the fasta file and rename to gene_name
        elif ref_name == line[1:]:
          name = gene_name
          seq = ''
          
        else: # get name of the test sequences   
          name = line[1:]
          seq = ''
        
      else: # sequence substring
        line = ''.join(line.split())
        seq += line.upper()
        # this is where the dictionary is defined
        if name == gene_name:
          if name in seq_dict.keys():
            # delete the dictionary entry with incomplete sequence
            # this happens when the ref sequence is in multiple lines
            del seq_dict[name]
          # ensure that the ref gene name and its sequence are at the top of the dictionary   
          temp_dict = {name:seq}
          # this will override dictionary items with the same key/name
          temp_dict.update(seq_dict)
          seq_dict = temp_dict

        else:
          seq_dict[name] = seq
  
  return seq_dict, ref_name

def create_sequences(seq_dict, start_loc, stop_loc):
  for name in seq_dict:
      seq = seq_dict[name]
      seq_dict[name] = Sequence(name, seq, start_loc, stop_loc) # instantiate Sequence

def calculate_protein_diffs(seq_dict, gene_name, ref_protein):
  for name in seq_dict:
    # trigger protein diff calculation if sequence is not reference
    seq = seq_dict[name]
    seq.calculate_protein_diffs(gene_name, ref_protein, name==gene_name)

# Not necessary
def parse_epitopes(epitope_file, seq_dict, ref_protein):
  with open(epitope_file) as eps:
    for line in eps:
      line = line.rstrip()
      values = line.split('-')
      start = int(values[0])
      end = int(values[1])

      for seq in seq_dict.values():
        # add epitope to sequence
        # calculate epitope protein diff against reference protein
        seq.add_epitope(start, end, ref_protein)

def search_start_stop(ref_genome, gene_start, gene_stop):
  start_loc = None
  stop_loc = None

  for i in range(1,len(ref_genome)):
    if ref_genome[i-1] == '-' and ref_genome[i] != '-':
      if ref_genome[i + gene_start : i + gene_start + 3] == 'ATG' and start_loc == None:
        # location of the first nucleotide
        start_loc = i + gene_start
        # location of the last nucleotide ['TAG','TGA','TAA']
        stop_loc = i + gene_stop
        return (start_loc, stop_loc)
        
    elif i == len(ref_genome) and (start_loc == None or stop_loc == None):
      # Raise Error Message
      print("Start and/or stop codon(s) not found in the reference sequence")

  return gene_start, gene_stop
      
  
def align_seq(inputs_fasta, reference_fasta):
  # Input the location of your mucle alignment program
  # Current available programs:
  # (1) muscle3.8.31_i86darwin64
  # (2) mafft [default]
  # MAFFT has the capability to align multiple sequences against a reference sequence https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
  muscle_exe = 'programs/muscle3.8.31_i86darwin64' 
  mafft_exe = 'programs/mafft'
  
  align_using = 'mafft'
  
  print('Aligning sequences...')
  if align_using == 'muscle':
    os.system(muscle_exe + ' -in ' + inputs_fasta + ' -out input/temp_aligned.fasta')
  else:
    os.system(mafft_exe + ' --quiet --keeplength --6merpair --addfragments ' + inputs_fasta + ' ' + reference_fasta + ' > input/temp_aligned.fasta' )
    
  print('Aligning sequences...[Completed]\n')
    
# Replace a character on a specific index of a string    
def replace_at(string, index, rep):    
  return string[0:index] + rep + string[index+1:]  
    
# Remove entire column of gaps according to the reference genome    
def remove_gaps(seq_dict, ref_name):
  # Search for location of the gaps in the reference genome
  indexes = [i for i, c in enumerate (seq_dict[ref_name]) if c == '-']
  
  # Delete the entire column of gap
  for key in seq_dict.keys():
    temp = seq_dict[key]
    for index in indexes:
      temp = replace_at(temp, index, '|')
      seq_dict[key] = temp.replace('|','')
  
  # Return the dictionary with deleted columns of gaps    
  return seq_dict

def parse_input_txt(input_txt):
  parsed = []
  with open(input_txt) as data:
    for line in data:
      line = line.rstrip() # remove trailing whitespaces
      line = line.split('-')
      parsed.append(line)
  return parsed

#######################################################################

if __name__ == "__main__":  
  
  # Welcome Message
  width = os.get_terminal_size().columns
  os.system('clear')
  print('\n\n\n')
  print('================================================================\n'.center(width))
  print('SARS-COV-2 EPITOPE SURVEILLANCE'.center(width))
  print('A Bioinformatics Project of Center for Informatics'.center(width))
  print('University of San Agustin, Iloilo City, Philippines'.center(width))
  print('Copyright \xa9 2021\n'.center(width))
  print('Author: Jonathan Adam A. Rico, MSc.\n'.center(width))
  print('================================================================\n\n'.center(width))

  
  # Align sequences using MUSCLE/MAFFT if the input is not yet aligned 
  if args.needs_alignment:
    align_seq(args.input, args.ref_genome)
    args.input = 'input/temp_aligned.fasta'
  # get reference genome name
  ref_genome_name = args.ref_name
  
  # Empty the output files
  open('output/protein.fasta', 'w').close() 
  open(args.out, 'w').close()
  
  # Create a list of B-cell epitope locations
  epitopes = parse_input_txt(args.eps)      
  # Create a list of gene names
  genes = parse_input_txt(args.gene_loc)
      
  for i in range(0,len(genes)):
    gene_name = genes[i][0]
    gene_start = int(genes[i][1]) - 1
    gene_stop = int(genes[i][2])
    
    # read and process file
    # create sequence dictionary
    seq_dict, ref_genome_name = parse_alignment(args.input, ref_genome_name, gene_name)

    # remove gaps in the reference genome and the corresponding column
    #seq_dict = remove_gaps(seq_dict, gene_name)

    # search for the locations of the ref genome sequence start and stop codon
    # since the genome is not yet trimmed, gene_name still contains the entire genome
    ref_genome = seq_dict[gene_name]
    (start_loc, stop_loc) = search_start_stop(ref_genome, gene_start, gene_stop)
    print("{} | start codon: {} | stop codon: {}\n".format(gene_name, start_loc+1, stop_loc))

    # Create a string of values for bcell epitope
    bcell_epitope = ''
    last_epi_stop = 0
    for j in range(len(epitopes)):
      if gene_name == epitopes[j][0]:
        epi_start = int(epitopes[j][1]) - 1
        epi_stop = int(epitopes[j][2])
        bcell_epitope = str(bcell_epitope) + '_' * (epi_start - last_epi_stop) + 'W' * (epi_stop - epi_start)
        #print('B-cell epitope:{} '.format(epi_stop - epi_start))
        last_epi_stop = epi_stop

    # Include the reference genome at the top of the Trimmed output file for checking purposes
    if i==0:
      # Initialize output files
      open('output/trimmed.fasta', 'w').write('>' + str(ref_genome_name) + '\n' + str(ref_genome) + '\n')  
    
    open('output/protein.fasta', 'a').write('> B-Cell Epitope Regions\n' + str(bcell_epitope) + '\n')
    open(args.out, 'a').write('> B-Cell Epitope Regions\n' + str(bcell_epitope) + '\n')

    # initialize sequences
    # seq_dict = { "name1": "SEQUENCE", "name2": "SEQUENCE" }
    create_sequences(seq_dict, start_loc, stop_loc)
    # seq_dict = { "name1": Sequence(), "name2": Sequence() }
  
    # calculate protein diffs for each test sequence against the reference protein
    ref_protein = seq_dict[gene_name].protein 
    calculate_protein_diffs(seq_dict, gene_name, ref_protein)

    # parse epitopes (Not necessary)
    # calculate epitope protein diffs for each sequence against reference protein
    #if args.eps is not None:
    #  parse_epitopes(args.eps, seq_dict, ref_protein)
      
    # Append results in a single output file
    with open(args.out, 'a') as out:
      for seq in seq_dict.values():
        # Write the reference sequence and test sequences to the output file
        if seq.name == gene_name or seq.name not in genes:
          out.write('>' + str(seq.name) + '\n') 
          out.write(str(''.join(seq.mutations)) + '\n')

          
