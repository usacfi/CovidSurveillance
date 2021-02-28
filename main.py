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
python3 main.py -o output/mutations.fasta -g input/Reference_genes_locations.txt -i input/010821-011521_indonesia.fasta
-ref "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome" -aln True
'''


parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-i', '--input_alignments_fasta', dest='input', help='input alignment file', required=True)

parser.add_argument('-aln', '--needs_alignment', dest='needs_alignment', help='(Type: Bool) True if input still needs to be aligned,'
                    'False if input is already aligned, e.g. -aln True, default value is False', 
                    default=False, type=lambda x: (str(x).lower() in ['true', '1', 'yes']))

parser.add_argument('-g', '--genes', dest='genes', help='input file with the location of genes', required=True)

parser.add_argument('-o', '--output', dest='out', help='base name for output files', required=True)

parser.add_argument('-ref', '--reference_name', dest='ref',
                    help='name of reference sequence in alignment to compare with other sequences; '
                         'otherwise uses first sequence in alignment as reference')

parser.add_argument('-eps', '--epitope_regions', dest='eps', help='text file that lists epitopes to check')

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
      
def plot_mutations(seq_dict, gene_name):
  
  arr = [v for v in seq_dict.values()]
  arr_flat = list(chain.from_iterable(arr))
  
  print(arr)
  print(arr_flat)
  
  amino = {'_' : 0, '-' : 1, 'F' : 2, 'C' : 3, 'A' : 4, 'I' : 5, 'E' : 6, 'D' : 7,
           'G' : 8, 'L' : 9, 'K' : 10, 'Q' : 11, 'W' : 12, 'Y' : 13, 'H' : 14, 'R' : 15,
           'T' : 16, 'N' : 17, 'M' : 18, 'V' : 19, 'P' : 20, 'S' : 21, 'X' : 22}
           
  df = pd.DataFrame(np.array([amino[i] for i in arr_flat]).reshape(15,len(arr)))    
  fig, ax = plt.subplots()   
  ax.imshow(df.values, vmin=0, cmap='jet')  
  
  for i in range(len(arr)):
    for j in range(len(arr[0])):
      ax.text(j,i, arr[i][j], ha='center', va='center')
  
  plt.title(str(gene_name))
  ax.set_yticklabels([])
  ax.set_xticklabels([])
  ax.set_yticks([])
  ax.set_xticks([])
  
def align_seq(fasta_file):
  # Input the location of your mucle alignment program
  # Current available programs:
  # (1) muscle3.8.31_i86darwin64
  # (2) mafft --default
  muscle_exe = 'programs/muscle3.8.31_i86darwin64' 
  mafft_exe = 'programs/mafft'
  
  align_using = 'mafft'
  
  if align_using == 'muscle':
    os.system(muscle_exe + ' -in ' + fasta_file + ' -out input/temp_aligned.fasta')
  else:
    os.system(mafft_exe + ' ' + fasta_file + ' > input/temp_aligned.fasta' )

#######################################################################

if __name__ == "__main__":
 
# Align sequences using MUSCLE if the input is not yet aligned 
  if args.needs_alignment:
    align_seq(args.input)
    args.input = 'input/temp_aligned.fasta'
  
  # Empty the output files
  init_outfiles = True
#  open('output/trimmed.fasta', 'w').close()
  open('output/protein.fasta', 'w').close() 
  open(args.out, 'w').close()
  
  # Create a list of gene names
  genes = []
  with open(args.genes) as gene_data:
    for line in gene_data:
      line = line.rstrip() # remove trailing whitespaces
      line = line.split('-')
      genes.append(line)
      
  for i in range(0,len(genes)):
    gene_name = genes[i][0]
    gene_start = int(genes[i][1]) - 1
    gene_stop = int(genes[i][2])

    # get reference genome name
    ref_genome_name = args.ref
    # read and process file
    # create sequence dictionary
    seq_dict, ref_genome_name = parse_alignment(args.input, ref_genome_name, gene_name)

    # search for the locations of the ref genome sequence start and stop codon
    # since the genome is not yet trimmed, gene_name still contains the entire genome
    ref_genome = seq_dict[gene_name]
    (start_loc, stop_loc) = search_start_stop(ref_genome, gene_start, gene_stop)
    print("start codon at {} \nstop codon at {}".format(start_loc, stop_loc))

    # Include the reference genome at the top of the Trimmed output file for checking purposes
    if init_outfiles:
      # Initialize output files
      open('output/trimmed.fasta', 'w').write('>' + str(ref_genome_name) + '\n' + str(ref_genome) + '\n')
      init_outfiles = False   

    # initialize sequences
    # seq_dict = { "name1": "SEQUENCE", "name2": "SEQUENCE" }
    create_sequences(seq_dict, start_loc, stop_loc)
    # seq_dict = { "name1": Sequence(), "name2": Sequence() }
  
    # calculate protein diffs for each test sequence against the reference protein
    ref_protein = seq_dict[gene_name].protein 
    calculate_protein_diffs(seq_dict, gene_name, ref_protein)

    # parse epitopes
    # calculate epitope protein diffs for each sequence against reference protein
    if args.eps is not None:
      parse_epitopes(args.eps, seq_dict, ref_protein)
      
    # Append results in a single output file
    with open(args.out, 'a') as out:
      for seq in seq_dict.values():
        # Write the reference sequence and test sequences to the output file
        if seq.name == gene_name or seq.name not in genes:
          out.write('>' + str(seq.name) + '\n') 
          out.write(str(''.join(seq.mutations)) + '\n')

          
    # Plot the protein differences
    #plot_mutations(seq_dict, gene_name)
          
  
  
  # TODO: output separate file with epitope information; different format?
