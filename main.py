# -*- coding:utf-8 -*-

import argparse
from classes.Sequence import Sequence
from classes.Functions import *
import pandas as pd



'''
Sample Run in Command Prompt/Terminal: 
python3 main.py -i references/Sequences/Morocco/48_Morocco_gisaid_hcov-19_2020_07_21_03.fasta -aln True
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
  
  # Create a dataframe with mutation and respective instances
  new_df = create_df(df, genes, epitopes)
  print(new_df[new_df.columns[0:6]])
  
  # Save the dataframes into separate files
  df.to_csv('output/05_aminoacid_replacements.csv')
  new_df.to_csv('output/06_mutation_profile.csv')
  
  # Plot the mutation profile
  plot_mutation_profile(new_df, genes, len(df))  
  
  # Plot mutations geographically
  geoplot_mutations(new_df)   
