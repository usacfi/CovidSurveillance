# -*- coding:utf-8 -*-

import argparse
from classes.Sequence import Sequence
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from itertools import chain
import seaborn as sns
import os

'''
Sample Run in Command Prompt/Terminal: 
python3 main.py -i references/Sequences/Morocco/48_Morocco_gisaid_hcov-19_2020_07_21_03.fasta -loc input/reference_genes_locations.txt 
-o output/04_mutations.fasta -ref "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome" -aln True
'''


parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-i', '--input_alignments_fasta', dest='input', help='input alignment file', required=True)

parser.add_argument('-aln', '--needs_alignment', dest='needs_alignment', help='(Type: Bool) True if input still needs to be aligned,'
                    'False if input is already aligned, e.g. -aln True, default value is True', 
                    default=True, type=lambda x: (str(x).lower() in ['true', '1', 'yes']))

parser.add_argument('-loc', '--gene_loc', dest='gene_loc', help='input file with the location of genes', default='input/reference_genes_locations.txt' )

parser.add_argument('-o', '--output', dest='out', help='base name for output files', default='output/04_mutations.fasta')

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

def calculate_protein_diffs(seq_dict, gene_name, ref_gene, start_loc, df):
  for name in seq_dict:
    # trigger protein diff calculation if sequence is not reference
    seq = seq_dict[name]
    df = seq.calculate_protein_diffs(name, gene_name, ref_gene, name==gene_name, start_loc, df)
  return df


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
    os.system(muscle_exe + ' -in ' + inputs_fasta + ' -out output/01_temp_aligned.fasta')
  else:
    os.system(mafft_exe + ' --quiet --keeplength --6merpair --addfragments ' + inputs_fasta + ' ' + reference_fasta + ' > output/01_temp_aligned.fasta' )
    
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

def is_radical(mutation):
  replacement = mutation[0] + mutation[-1]
  reverse = mutation[-1] + mutation[0]
  
  # Based on Grantham's mean chemical difference index (1974)
  criteria_table = ['RS', 'LS', 'LR', 'PR', 'AR', 'VR', 'GR', 'GL', 'GV', 'IS', 'IG', 'FS',
                    'FP', 'FT', 'FA', 'FG', 'YS', 'YP', 'YA', 'YG', 'CS', 'CR', 'CL', 'CP',
                    'CT', 'CA', 'CV', 'CG', 'CI', 'CF', 'CY', 'HC', 'QL', 'QI', 'QF', 'QC',
                    'NL', 'NA', 'NV', 'NI', 'NF', 'NY', 'NC', 'KS', 'KL', 'KP', 'KA', 'KG',
                    'KI', 'KF', 'KC', 'DL', 'DP', 'DA', 'DV', 'DI', 'DF', 'DY', 'DC', 'DK',
                    'EL', 'EA', 'EV', 'EI', 'EF', 'EY', 'EC', 'MS', 'MG', 'MC', 'MQ', 'MN',
                    'MD', 'ME', 'WS', 'WR', 'WP', 'WT', 'WA', 'WG', 'WC', 'WH', 'WQ', 'WN', 
                    'WK', 'WD', 'WE']
  
  if replacement in criteria_table:
    return 'Radical'
  elif reverse in criteria_table:
    return 'Radical'
  else:
    return 'Conservative'
    
def instances_per_region(df, site, mutation):
  
  instance_list = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]              
  location = df[df[site]==mutation].Location

  # Enumerate all possible names of the regions
  region_dict = { 0  : ['NCR', 'Manila', 'National Capital Region'], 
                  1  : ['CAR', 'Cordillera Administrative Region'],
                  2  : ['Region I', 'Ilocos Region'],
                  3  : ['Region II','Cagayan Valley'],
                  4  : ['Region III', 'Central Luzon'],
                  5  : ['Region IV-A', 'Calabarzon'],
                  6  : ['Region IV-B', 'Mimaropa'],
                  7  : ['Region V', 'Bicol Region'],
                  8  : ['Region VI', 'Western Visayas'],
                  9  : ['Region VII', 'Central Visayas'],
                  10 : ['Region VIII', 'Eastern Visayas'],
                  11 : ['Region IX','Zamboanga Peninsula'],
                  12 : ['Region X','Northern Mindanao'],
                  13 : ['Region XI','Davao Region'],
                  14 : ['Region XII', 'SOCCKSARGEN'],
                  15 : ['Caraga', 'Davao Oriental'],
                  16 : ['BARMM', 'ARMM']}
           
  for key in region_dict.keys():
    for i in location.index.values:
      loc_name = location[i].replace('Asia / Philippines / ', '')
      if loc_name in region_dict[key]:
        instance_list[key] += 1
  
  return instance_list          
  
  
def plot_mutation_profile(df, genes, len_inputs):
  g = sns.barplot(x=df['mutation'], y=df['instances'], hue=df['gene'], saturation=0.5, dodge=False)
  g.legend(loc='upper center', ncol=len(genes), title='Gene')

  # Annotate if the mutations fall under epitope regions
  for i in range(len(df['mutation'])):
    if df['replacement'][i] == 'Radical':
      color = 'red'
    else:
      color = 'green'
    
    if df['epitope'][i] != None:
      g.text(i, df['instances'][i]+0.25, df['epitope'][i][0], fontdict=dict(color=color, fontsize=8),
             horizontalalignment='center')
  plt.title('Mutation profile of {} Philippine Sars-Cov-2 samples'.format(len_inputs-1))           
  plt.ylabel('Total # of mutation instances (out of {})'.format(len_inputs-1))
  plt.xticks(rotation=60, fontsize=6, ha='right')
  plt.tight_layout()
  plt.savefig('output/07_mutation_profile.pdf', orientation='landscape', format='pdf')

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
  print('Authors: Rico JA, Zamora PRF, Bolinas DK, Aguila RN,'.center(width))
  print('Isip I, Dumancas G, Fenger D, de Castro R\n'.center(width))
  print('================================================================\n\n'.center(width))

  # Open Metadata file
  meta_df = pd.read_csv('references/Sequences/Philippines/Other Files/[Patient Status Metadata] gisaid_hcov-19_2021_03_23_08.tsv', sep='\t')
  # Combine three columns into one column to match the ID
  meta_df['ID'] = meta_df[meta_df.columns[0:3]].apply(
    lambda x: '|'.join(x.dropna().astype(str)), axis=1
  )
  
  # Align sequences using MUSCLE/MAFFT if the input is not yet aligned 
  if args.needs_alignment:
    align_seq(args.input, args.ref_genome)
    args.input = 'output/01_temp_aligned.fasta'
  # get reference genome name
  ref_genome_name = args.ref_name
  
  # Empty the output files
  open('output/03_protein.fasta', 'w').close() 
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

    # Create a string of values for each epitopes
    bcell_epitope = '_' * int((gene_stop - gene_start + 1)/3 - 1) 
    tcell_epitope = '_' * int((gene_stop - gene_start + 1)/3 - 1) 
    for j in range(len(epitopes)):
      if epitopes[j][0] == 'B cell':
        if gene_name == epitopes[j][1]:
          epi_start = int(epitopes[j][2]) - 1
          epi_stop = int(epitopes[j][3])
          bcell_epitope = bcell_epitope[0:epi_start] + "W" * (epi_stop - epi_start) + bcell_epitope[epi_stop:]
          #print('B-cell epitope:{} '.format(epi_stop - epi_start))
      elif epitopes[j][0] == 'T cell':
        if gene_name == epitopes[j][1]:
          epi_start = int(epitopes[j][2]) - 1
          epi_stop = int(epitopes[j][3])
          tcell_epitope = tcell_epitope[0:epi_start] + "W" * (epi_stop - epi_start) + tcell_epitope[epi_stop:]
          #print('T-cell epitope:{} {} '.format(epi_start+1, epi_stop))          
        

    # Include the reference genome at the top of the Trimmed output file for checking purposes
    if i==0:
      # Initialize output files
      open('output/02_trimmed.fasta', 'w').write('>' + str(ref_genome_name) + '\n' + str(ref_genome) + '\n')
      df = pd.DataFrame(seq_dict.keys(), columns=['ID'])  
      df['ID'] = df['ID'].replace([gene_name],ref_genome_name)
      df = df.merge(meta_df, how='left', on='ID')
    
    open('output/03_protein.fasta', 'a').write('> B-Cell Epitope Regions\n' + str(bcell_epitope) + '\n')
    open('output/03_protein.fasta', 'a').write('> T-Cell Epitope Regions\n' + str(tcell_epitope) + '\n')    
    open(args.out, 'a').write('> B-Cell Epitope Regions\n' + str(bcell_epitope) + '\n')
    open(args.out, 'a').write('> T-Cell Epitope Regions\n' + str(tcell_epitope) + '\n')    

    # initialize sequences
    # seq_dict = { "name1": "SEQUENCE", "name2": "SEQUENCE" }
    create_sequences(seq_dict, start_loc, stop_loc)
    # seq_dict = { "name1": Sequence(), "name2": Sequence() }
  
    # calculate protein diffs for each test sequence against the reference protein
    ref_gene = seq_dict[gene_name].protein 
    df = calculate_protein_diffs(seq_dict, gene_name, ref_gene, start_loc, df)

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
  
  # Visualization of Variants vs. Number of Occurrences
  aa_diffs = []        
  for site in df.columns.values:
    if isinstance(site, int):
      for mutation in df[site].dropna().unique():
        radical_replacement = is_radical(mutation)
        region_instances = instances_per_region(df, site, mutation) 
        # Create gene and epitope columns
        for i in range(len(genes)):
          if site >= int(genes[i][1])-1 and site <= int(genes[i][2]):
            gene = genes[i][0]
            epitope = None
            for j in range(len(epitopes)):
              if gene == epitopes[j][1]:
                if site >= int(genes[i][1])-1 + int(epitopes[j][2])-1 and site <= int(genes[i][1])-1 + int(epitopes[j][3]):
                  epitope = epitopes[j][0]

            
        count = len(df[df[site]==mutation])
        cols = [site, mutation, count, gene, epitope, radical_replacement]
        cols.extend(region_instances)
        aa_diffs.append(cols)
                
  new_df = pd.DataFrame(aa_diffs, columns=['site', 'mutation','instances', 'gene', 'epitope', 'replacement',
                                            'NCR', 'CAR', 'Region I', 'Region II', 'Region III', 'Region IV-A',
                                            'Region IV-B', 'Region V', 'Region VI', 'Region VII', 'Region VIII',
                                            'Region IX', 'Region X', 'Region XI', 'Region XII', 'Caraga', 'BARMM'])
  new_df = new_df.sort_values(by='site')
  new_df = new_df.reset_index(drop=True)
  print(new_df[new_df.columns[0:6]])
  
  # Save the dataframes into separate files
  df.to_csv('output/05_aminoacid_replacements.csv')
  new_df.to_csv('output/06_mutation_profile.csv')
  
  # Plot the mutation profile
  plot_mutation_profile(new_df, genes, len(df))     
