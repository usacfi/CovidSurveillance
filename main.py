import argparse
from classes.Sequence import Sequence
import numpy as np

parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-aln', '--alignments_fasta', dest='aln', help='input alignment file', required=True)

parser.add_argument('-o', '--output', dest='out', help='base name for output files', required=True)

parser.add_argument('-ref', '--reference_name', dest='ref',
                    help='name of reference sequence in alignment to compare with other sequences; '
                         'otherwise uses first sequence in alignment as reference')

parser.add_argument('-eps', '--epitope_regions', dest='eps', help='text file that lists epitopes to check')

args = parser.parse_args()

def parse_alignment(fasta_file):
  ref_name = args.ref
  seq_dict = dict()

  with open(fasta_file) as data:
    name = ''
    seq = ''

    for line in data:
      line = line.rstrip() # remove trailing whitespaces

      if line.startswith('>'): # name
        if ref_name == None: # if user did not specify -ref input
          ref_name = line[1:] # assumes first name found is the reference
          
        name = line[1:]
        seq = ''
        
      else: # sequence substring
        line = ''.join(line.split())
        seq += line.upper()
        # this is where the dictionary is defined
        seq_dict[name] = seq

  
  return (ref_name, seq_dict)

def create_sequences(seq_dict):
  for name in seq_dict:
    seq = seq_dict[name]
    seq_dict[name] = Sequence(name, seq) # instantiate Sequence

def calculate_protein_diffs(seq_dict, ref_name, ref_protein):
  for name in seq_dict:
    if name != ref_name:
      # trigger protein diff calculation if sequence is not reference
      seq = seq_dict[name]
      seq.calculate_protein_diffs(ref_name, ref_protein)

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

def search_start_stop(ref_protein):
  start_loc = None
  stop_loc = None
  
  for i in range(len(ref_protein)):
    if ref_protein[i] == 'M' and start_loc == None:
      start_loc = i
    elif ref_protein[i] == '-' and start_loc != None and stop_loc == None:
      stop_loc = i-1
      return start_loc, stop_loc
    elif i == len(ref_protein) and (start_loc == None or stop_loc == None):
      # Raise Error Message
      print("Start and stop codons not found in the reference sequence")
    

#######################################################################

if __name__ == "__main__":
  # read and process file
  # get reference name
  # create sequence dictionary
  (ref_name, seq_dict) = parse_alignment(args.aln)

  print("reference name: " + ref_name)


  # initialize sequences
  # seq_dict = { "name1": "SEQUENCE", "name2": "SEQUENCE" }
  create_sequences(seq_dict)
  # seq_dict = { "name1": Sequence(), "name2": Sequence() }
  
  # define the reference protein sequence
  ref_protein = seq_dict[ref_name].protein

  # search for the locations of the ref seq's start and stop codon
  start_loc, stop_loc = search_start_stop(ref_protein)
  #print("start codon at {} \nstop codon at {}".format(start_loc, stop_loc))

  # calculate protein diffs for each sequence against the reference protein
  calculate_protein_diffs(seq_dict, ref_name, ref_protein)

  # parse epitopes
  # calculate epitope protein diffs for each sequence against reference protein
  if args.eps is not None:
    parse_epitopes(args.eps, seq_dict, ref_protein)

  with open(args.out + '.xls', 'w') as out:
    for seq in seq_dict.values():
      out.write(str(seq) + '\n')
      #print(str(seq))

  # TODO: output separate file with epitope information; different format?
