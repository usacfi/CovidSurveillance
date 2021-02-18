import argparse
from classes.Sequence import Sequence
import numpy as np

'''
Sample Run in Command Prompt/Terminal: 
python3 main.py -aln input/010821-011521_indonesia.fasta -o output/test_output.xls -genes input/Reference_genes_locations.txt
'''


parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-aln', '--alignments_fasta', dest='aln', help='input alignment file', required=True)

parser.add_argument('-genes', '--genes', dest='genes', help='input file with the location of genes', required=True)

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
          
        else:    
          name = line[1:]
          seq = ''
        
      else: # sequence substring
        line = ''.join(line.split())
        seq += line.upper()
        # this is where the dictionary is defined
        seq_dict[name] = seq
  
  return seq_dict

def create_sequences(seq_dict, start_loc, stop_loc):
  for name in seq_dict:
    seq = seq_dict[name]
    seq_dict[name] = Sequence(name, seq, start_loc, stop_loc) # instantiate Sequence

def calculate_protein_diffs(seq_dict, ref_name, ref_protein):
  for name in seq_dict:
    # trigger protein diff calculation if sequence is not reference
    seq = seq_dict[name]
    seq.calculate_protein_diffs(ref_name, ref_protein, name==ref_name)

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

def search_start_stop(ref_seq, gene_start, gene_stop):
  start_loc = None
  stop_loc = None

  for i in range(1,len(ref_seq)):
    if ref_seq[i-1] == '-' and ref_seq[i] != '-':
      if ref_seq[i + gene_start : i + gene_start + 3] == 'ATG' and start_loc == None:
        # location of the first nucleotide
        start_loc = i + gene_start
        # location of the last nucleotide ['TAG','TGA','TAA']
        stop_loc = i + gene_stop
        return (start_loc, stop_loc)
        
    elif i == len(ref_seq) and (start_loc == None or stop_loc == None):
      # Raise Error Message
      print("Start and/or stop codon(s) not found in the reference sequence")

  return gene_start, gene_stop
      

#######################################################################

if __name__ == "__main__":
  
  # Empty the output file
  open(args.out, 'w').close()
  
  # Create a list of reference names
  genes = []
  with open(args.genes) as gene_data:
    for line in gene_data:
      line = line.rstrip() # remove trailing whitespaces
      line = line.split('-')
      genes.append(line)
      
  
  ref_name = args.ref
  for i in range(0,len(genes)):
    gene_name = genes[i][0]
    gene_start = int(genes[i][1]) - 1
    gene_stop = int(genes[i][2])

    # read and process file
    # get reference name
    # create sequence dictionary
    seq_dict = parse_alignment(args.aln, ref_name, gene_name)

    # search for the locations of the ref seq's start and stop codon
    ref_seq = seq_dict[gene_name]
    (start_loc, stop_loc) = search_start_stop(ref_seq, gene_start, gene_stop)
    print("start codon at {} \nstop codon at {}".format(start_loc, stop_loc))

    # initialize sequences
    # seq_dict = { "name1": "SEQUENCE", "name2": "SEQUENCE" }
    create_sequences(seq_dict, start_loc, stop_loc)
    # seq_dict = { "name1": Sequence(), "name2": Sequence() }
  
    # calculate protein diffs for each sequence against the reference protein
    ref_protein = seq_dict[gene_name].protein 
    calculate_protein_diffs(seq_dict, gene_name, ref_protein)

    # parse epitopes
    # calculate epitope protein diffs for each sequence against reference protein
    if args.eps is not None:
      parse_epitopes(args.eps, seq_dict, ref_protein)

    # Append results in a single output file
    with open(args.out, 'a') as out:
      for seq in seq_dict.values():
          # Write the reference sequence to the output file
          if seq.name == gene_name:
            out.write(str(seq) + '\n')
            #print(str(seq))
          # Write the input sequences to the output file
          elif seq.name not in genes:
            out.write(str(seq) + '\n')
            #print(str(seq))
          

  # TODO: output separate file with epitope information; different format?
