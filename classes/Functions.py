from classes.Sequence import Sequence
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import seaborn as sns
import os
import osmnx as ox
import time
import copy



def parse_alignment(fasta_file, ref_name, gene_name):
  ''' Creates a dataframe from the fasta file and also returns the reference name '''
  
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



def fasta_to_df(fasta_file):
  ''' Converts a fasta file into a data frame '''
  
  with open(fasta_file) as data:
    name = []
    sequence = []
    for line in data:
      line = line.rstrip()
      if line.startswith('>'):
        name.append(line[1:])
      else:
        sequence.append(line)
        
  return pd.DataFrame({'name':name,'sequence':sequence})      



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
      
  
  
def align_seq(inputs_fasta, reference_fasta, output_fasta):
  # Input the location of your mucle alignment program
  # Current available programs:
  # (1) muscle3.8.31_i86darwin64
  # (2) mafft [default]
  # MAFFT has the capability to align multiple sequences against a reference sequence https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
  muscle_exe = 'programs/muscle3.8.31_i86darwin64' 
  mafft_exe = 'programs/mafft'
  
  align_using = 'mafft'
  
  print('Aligning sequences...')
  tic = time.perf_counter()
  if align_using == 'muscle':
    os.system(muscle_exe + ' -in ' + inputs_fasta + ' -out ' + output_fasta)
  else:
    os.system(mafft_exe + ' --quiet --keeplength --6merpair --addfragments ' + inputs_fasta + ' ' + reference_fasta + ' > ' + output_fasta )
  toc = time.perf_counter()  
  print('Aligning sequences...[Completed]')
  print('Alignment Runtime: {:.4f} seconds\n'.format(toc-tic))
    
    
     
def replace_at(string, index, rep):   
  ''' Replace a character on a specific index of a string ''' 
  return string[0:index] + rep + string[index+1:]  
    

def replace_with(my_lists, old, new):
  ''' Replace a character in a list of lists or numpy array'''
  for i in range(len(my_lists)):
    my_lists[i] = [new if x==old else x for x in my_lists[i]]
    
  return my_lists
    
    
# Not necessary 
def remove_gaps(seq_dict, ref_name):
  ''' Delete the entire column of gaps according to the reference genome '''
  
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
  ''' Converts a text file into a list of lists
    The text file (.txt) values should be separated by a dash "-" and new line for another list '''
  
  parsed = []
  with open(input_txt) as data:
    for line in data:
      line = line.rstrip() # remove trailing whitespaces
      line = line.split('-')
      parsed.append(line)
  return parsed



def is_radical(mutation):
  ''' Classifies a mutation into 1 (radical) or 0 (conservative) '''
  
  replacement = mutation[0] + mutation[-1]
  reverse = mutation[-1] + mutation[0]
  
  # Based on Grantham's mean chemical difference index (1974)
  grantham = ['RS', 'LS', 'LR', 'PR', 'AR', 'VR', 'GR', 'GL', 'GV', 'IS', 'IG', 'FS',
              'FP', 'FT', 'FA', 'FG', 'YS', 'YP', 'YA', 'YG', 'CS', 'CR', 'CL', 'CP',
              'CT', 'CA', 'CV', 'CG', 'CI', 'CF', 'CY', 'HC', 'QL', 'QI', 'QF', 'QC',
              'NL', 'NA', 'NV', 'NI', 'NF', 'NY', 'NC', 'KS', 'KL', 'KP', 'KA', 'KG',
              'KI', 'KF', 'KC', 'DL', 'DP', 'DA', 'DV', 'DI', 'DF', 'DY', 'DC', 'DK',
              'EL', 'EA', 'EV', 'EI', 'EF', 'EY', 'EC', 'MS', 'MG', 'MC', 'MQ', 'MN',
              'MD', 'ME', 'WS', 'WR', 'WP', 'WT', 'WA', 'WG', 'WC', 'WH', 'WQ', 'WN', 
              'WK', 'WD', 'WE']
  
  if replacement in grantham:
    return 'Radical'
  elif reverse in grantham:
    return 'Radical'
  else:
    return 'Conservative'


    
def radical_list(criteria):
  
  if criteria == 'grantham':
    radicals = ['AC', 'AD', 'AE', 'AF', 'AK', 'AN', 'AR', 'AW', 'AY', 'CA', 
                'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 
                'CP', 'CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY', 'DA', 'DC', 
                'DF', 'DI', 'DK', 'DL', 'DM', 'DP', 'DV', 'DW', 'DY', 'EA', 
                'EC', 'EF', 'EI', 'EL', 'EM', 'EV', 'EW', 'EY', 'FA', 'FC', 
                'FD', 'FE', 'FG', 'FK', 'FN', 'FP', 'FQ', 'FS', 'FT', 'GC', 
                'GF', 'GI', 'GK', 'GL', 'GM', 'GR', 'GV', 'GW', 'GY', 'HC', 
                'HW', 'IC', 'ID', 'IE', 'IG', 'IK', 'IN', 'IQ', 'IS', 'KA', 
                'KC', 'KD', 'KF', 'KG', 'KI', 'KL', 'KP', 'KS', 'KW', 'LC', 
                'LD', 'LE', 'LG', 'LK', 'LN', 'LQ', 'LR', 'LS', 'MC', 'MD', 
                'ME', 'MG', 'MN', 'MQ', 'MS', 'NA', 'NC', 'NF', 'NI', 'NL', 
                'NM', 'NV', 'NW', 'NY', 'PC', 'PD', 'PF', 'PK', 'PR', 'PW', 
                'PY', 'QC', 'QF', 'QI', 'QL', 'QM', 'QW', 'RA', 'RC', 'RG', 
                'RL', 'RP', 'RS', 'RW', 'SC', 'SF', 'SI', 'SK', 'SL', 'SM', 
                'SR', 'SV', 'SW', 'SY', 'TC', 'TF', 'TW', 'VC', 'VD', 'VE', 
                'VG', 'VN', 'VS', 'WA', 'WC', 'WD', 'WE', 'WG', 'WH', 'WK', 
                'WN', 'WP', 'WQ', 'WR', 'WS', 'WT', 'YA', 'YC', 'YD', 'YE', 
                'YG', 'YN', 'YP', 'YS']
                    
  elif criteria == 'blosum62':
    radicals = ['AD', 'AF', 'AH', 'AN', 'AW', 'AY', 'CD', 'CE', 'CF', 'CG', 
                'CH', 'CK', 'CN', 'CP', 'CQ', 'CR', 'CW', 'CY', 'DA', 'DC', 
                'DF', 'DI', 'DL', 'DM', 'DR', 'DV', 'DW', 'DY', 'EC', 'EF', 
                'EG', 'EI', 'EL', 'EM', 'EV', 'EW', 'EY', 'FA', 'FC', 'FD', 
                'FE', 'FG', 'FK', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'GC', 
                'GE', 'GF', 'GH', 'GI', 'GK', 'GL', 'GM', 'GP', 'GQ', 'GR', 
                'GT', 'GV', 'GW', 'GY', 'HA', 'HC', 'HG', 'HI', 'HL', 'HM', 
                'HP', 'HT', 'HV', 'HW', 'ID', 'IE', 'IG', 'IH', 'IK', 'IN', 
                'IP', 'IQ', 'IR', 'IS', 'IW', 'KC', 'KF', 'KG', 'KI', 'KL', 
                'KV', 'KW', 'KY', 'LD', 'LE', 'LG', 'LH', 'LK', 'LN', 'LP', 
                'LQ', 'LR', 'LS', 'LW', 'MD', 'ME', 'MG', 'MH', 'MN', 'MP', 
                'NA', 'NC', 'NF', 'NI', 'NL', 'NM', 'NP', 'NV', 'NW', 'NY', 
                'PC', 'PF', 'PG', 'PH', 'PI', 'PL', 'PM', 'PN', 'PR', 'PV', 
                'PW', 'PY', 'QC', 'QF', 'QG', 'QI', 'QL', 'QV', 'QW', 'RC', 
                'RD', 'RF', 'RG', 'RI', 'RL', 'RP', 'RV', 'RW', 'RY', 'SF', 
                'SI', 'SL', 'SV', 'SW', 'SY', 'TF', 'TG', 'TH', 'TW', 'TY', 
                'VD', 'VE', 'VG', 'VH', 'VK', 'VN', 'VP', 'VQ', 'VR', 'VS', 
                'VW', 'WA', 'WC', 'WD', 'WE', 'WG', 'WH', 'WI', 'WK', 'WL', 
                'WN', 'WP', 'WQ', 'WR', 'WS', 'WT', 'WV', 'XC', 'XP', 'XW', 
                'YA', 'YC', 'YD', 'YE', 'YG', 'YK', 'YN', 'YP', 'YR', 'YS', 
                'YT']
    
    ''' when blosum62 score is -3 and -4            
    radicals = ['AW', 'CD', 'CE', 'CG', 'CH', 'CK', 'CN', 'CP', 'CQ', 'CR', 
                'DC', 'DF', 'DI', 'DL', 'DM', 'DV', 'DW', 'DY', 'EC', 'EF', 
                'EI', 'EL', 'EW', 'FD', 'FE', 'FG', 'FK', 'FN', 'FP', 'FQ', 
                'FR', 'GC', 'GF', 'GI', 'GL', 'GM', 'GV', 'GY', 'HC', 'HI', 
                'HL', 'HV', 'ID', 'IE', 'IG', 'IH', 'IK', 'IN', 'IP', 'IQ', 
                'IR', 'IW', 'KC', 'KF', 'KI', 'KW', 'LD', 'LE', 'LG', 'LH', 
                'LN', 'LP', 'MD', 'MG', 'NC', 'NF', 'NI', 'NL', 'NV', 'NW', 
                'PC', 'PF', 'PI', 'PL', 'PW', 'PY', 'QC', 'QF', 'QI', 'RC', 
                'RF', 'RI', 'RV', 'RW', 'SW', 'VD', 'VG', 'VH', 'VN', 'VR', 
                'VW', 'WA', 'WD', 'WE', 'WI', 'WK', 'WN', 'WP', 'WR', 'WS', 
                'WV', 'YD', 'YG', 'YP']            
    '''
  else:
    print('Error: Specify your criteria on classfiying radical and conservative substitutions.'
          'Choose from: \n1. grantham \n2.blosum62 \n3.')
  
  return radicals
    
# Not necessary    
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
  

def create_df(df, genes, epitopes):
  aa_diffs = []        
  for site in df.columns.values:
    if isinstance(site, int):
      for mutation in df[site].dropna().unique():
        radical_replacement = is_radical(mutation)
        if 'Location' in df.columns: region_instances = instances_per_region(df, site, mutation) 
        count = len(df[df[site]==mutation])
        # Identify gene and epitope for the respective mutation
        for i in range(len(genes)):
          if site >= int(genes[i][1])-1 and site <= int(genes[i][2]):
            gene = genes[i][0]
            epitope = None
            for j in range(len(epitopes)):
              if gene == epitopes[j][1]:
                if site >= int(genes[i][1])-1 + int(epitopes[j][2])-1 and site <= int(genes[i][1])-1 + int(epitopes[j][3]):
                  epitope = epitopes[j][0]
      
        cols = [site, mutation, count, gene, epitope, radical_replacement]
        if 'Location' in df.columns: cols.extend(region_instances)
        aa_diffs.append(cols)
                
  ph_regions = ['NCR', 'CAR', 'Region I', 'Region II', 'Region III', 'Region IV-A',
                'Region IV-B', 'Region V', 'Region VI', 'Region VII', 'Region VIII',
                'Region IX', 'Region X', 'Region XI', 'Region XII', 'Caraga', 'BARMM']              
  column_names = ['site', 'mutation','instances', 'gene', 'epitope', 'replacement']
  
  if 'Location' in df.columns: column_names.extend(ph_regions)              
  new_df = pd.DataFrame(aa_diffs, columns=column_names)
  new_df = new_df.sort_values(by='site')
  new_df = new_df.reset_index(drop=True)
  
  return new_df

  
def which_protein(site, proteins):
  ''' proteins is a list of lists, e.g.
      [['Surface',123,456],
       ['Envelope',432,765],
       ['Membrane',234,567]] '''
  
  for protein in proteins:
    if site in range(int(protein[1]), int(protein[2])):
      return protein[0]   
  raise ValueError('Site location {} is outside the (nucleotide) regions of the relevant proteins'.format(site))     
  
  
def epi_hits(mutation, protein, epitopes):
  ''' Returns a list of epitopes hit given the mutation and protein '''
  
  epitopes_hit = []
  for epitope in epitopes:  
    if int(mutation[1:len(mutation)-1]) in range(int(epitope[2]), int(epitope[3])) and protein == epitope[1]:
      #epitopes_hit = epitopes_hit + epitope[0] + ','
      epitopes_hit.append(epitope[0])

  if len(epitopes_hit) == 0:
    epitopes_hit = 'None'
        
  return epitopes_hit  
  
  
# Not necessary  
def plot_mutation_profile(df, proteins, epitopes):
  tic = time.perf_counter()
  mutations = []
  mu_site = []
  # Get unique mutations per protein for the x-axis
  sites = list(df.iloc[:,1:].columns) # all columns
  for site in sites:
    prev_len = len(mutations)
    mutations.extend(df[site].unique())
    mu_site.extend([x for x in [site] * (len(mutations) - prev_len)])
  
  mu_df = pd.DataFrame({'mu_site':mu_site, 'mutations':mutations})  
  mu_df = mu_df.dropna()
  mu_df = mu_df.sort_values(by='mu_site')
  mu_df = mu_df.reset_index(drop=True)



  # Create a long stack of mutations 
  mu_stack = df[sites[0]]
  for site in sites[1:]:
    mu_stack.append(df[site])

  # Remove NaNs  
  mu_stack = mu_stack.dropna()  




  replacement_type = []
  mu_stack2int = []
  site_stack = []
  protein_stack = []
  epitope_stack = []
  for mu in list(mu_stack):
    # Convert into integers
    mu_stack2int.append(list(mu_df['mutations']).index(mu))
    # Classify mutations into radical and conservative
    replacement_type.append(is_radical(mu))
    # Site location of each mutation
    site = int(mu_df[mu_df['mutations']==mu].mu_site)
    site_stack.append(site)
    # Protein region of each mutation
    protein = which_protein(site, proteins)
    protein_stack.append(protein)
    # Epitopes hit
    epitope_stack.append(epi_hits(mu, protein, epitopes))


  stack_df = pd.DataFrame({'site':site_stack,
                          'mutation':mu_stack, 
                          'X_value':mu_stack2int, 
                          'Y_value':mu_stack.index, 
                          'replacement_type':replacement_type,
                          'protein':protein_stack,
                          'epitope':epitope_stack})
  stack_df = stack_df.sort_values(by='X_value')                        
  stack_df = stack_df.reset_index(drop=True)
  stack_df['epitope_str'] = stack_df.epitope.apply(', '.join)

  # Insert a 'radicality' column 
  stack_df.loc[stack_df.replacement_type == 'Radical', 'radicality'] = 1
  stack_df.loc[stack_df.replacement_type == 'Conservative', 'radicality'] = 0

  # Initialize arrays and lists with empty or NaN values
  radicality = np.full((len(df), len(mu_df)), np.nan)
  replacement_type = [x[:] for x in [[''] * len(mu_df)] * len(df)]
  mu_protein = [x[:] for x in [[''] * len(mu_df)] * len(df)]
  mu_epitope = [x[:] for x in [[''] * len(mu_df)] * len(df)]

  # Fill the arrays and lists
  for i in range(0,len(stack_df)):
    y = stack_df['Y_value'][i]
    x = stack_df['X_value'][i]
    radicality[y][x] = stack_df['radicality'][i]
    replacement_type[y][x] = stack_df['replacement_type'][i]
    mu_protein[y][x] = stack_df['protein'][i]
    mu_epitope[y][x] = stack_df['epitope'][i]  

  cnt_mu = []
  cnt_cons = []
  cnt_rad = []
  cnt_epi = []
  rad_score = []
  for k in df.index:
    count_mutations = stack_df.query('Y_value == {}'.format(k)).mutation.count()
    count_conservative = stack_df.query('Y_value == {} & replacement_type.str.contains("Conservative")'.format(k)).replacement_type.count()
    count_radical = stack_df.query('Y_value == {} & replacement_type.str.contains("Radical")'.format(k)).replacement_type.count()
    count_bcell = stack_df.query('Y_value == {} & epitope_str.str.contains("B cell")'.format(k)).epitope.count()
    count_tcell1 = stack_df.query('Y_value == {} & epitope_str.str.contains("T cell HLA class I")'.format(k)).epitope.count()
    count_tcell2 = stack_df.query('Y_value == {} & epitope_str.str.contains("T cell HLA class II")'.format(k)).epitope.count() 
    if count_mutations != 0:
      radicality_score = count_radical/count_mutations
    else:
      radicality_score = 0
    
    cnt_mu.append(count_mutations)
    cnt_cons.append(count_conservative)
    cnt_rad.append(count_radical)
    cnt_epi.append(count_bcell + count_tcell1 + count_tcell2)
    rad_score.append(round(radicality_score, 2))
  
  table_df = pd.DataFrame({ 'sample_id':df.index,
                            'count_mutation':cnt_mu, 
                            'count_conservative':cnt_cons, 
                            'count_radical':cnt_rad,
                            'count_epitope':cnt_epi,
                            'radicality_score':rad_score})


  table_plot = go.Figure(data=[go.Table(
      header=dict(values=list(table_df.columns),
                  fill_color='paleturquoise',
                  align='left'),
      cells=dict(values=[table_df.sample_id, 
                        table_df.count_mutation, 
                        table_df.count_conservative, 
                        table_df.count_radical, 
                        table_df.count_epitope, 
                        table_df.radicality_score],
                 fill_color='lavender',
                 align='left')) ])
                 
  table_plot.write_html('output/07.1_table.html')

  fig = go.Figure(data=go.Heatmap(
            x = mu_df['mutations'],
            y = list(df.index),
            z = radicality,
            customdata = np.dstack((replacement_type,mu_protein,mu_epitope)), 
            type = 'heatmap',
            colorscale = [[0,'green'],[1,'red']],
            xgap = 0.5,
            ygap = 0.5,
            hoverongaps=False,
            showscale=False,
            hovertemplate="Sample ID: %{y}<br>"
                          "Mutation: %{x}<br>"
                          "Type: %{customdata[0]}<br>"
                          "Protein: %{customdata[1]}<br>"
                          "Epitopes Hit: %{customdata[2]}<extra></extra>"))

  # Plot Protein Regions
  protein_colors = ["LightSkyBlue","Orange","RoyalBlue","Yellow"]
  for i in range(len(proteins)):
    start = list(stack_df[stack_df['protein']==proteins[i][0]].X_value)[0]
    end = list(stack_df[stack_df['protein']==proteins[i][0]].X_value)[-1]  
    fig.add_shape(type="rect",
        x0=start-0.5, y0=-5, x1=end+0.5, y1=-0.5,
        line=dict(width=0),fillcolor=protein_colors[i])
    fig.add_annotation(x=(start+end)/2, y=-3,text=proteins[i][0],showarrow=False) 
 

  fig.add_shape(type="rect",x0=-0.5, y0=len(df)+0.5, x1=len(mu_df)-0.5, y1=len(df)+10,
          line=dict(width=0),fillcolor="White")
  fig.add_shape(type="line",x0=-0.5, y0=len(df)+9.5, x1=len(mu_df)-0.5, y1=len(df)+9.5,line=dict(color='Grey'))
  fig.add_shape(type="line",x0=-0.5, y0=len(df)+6.5, x1=len(mu_df)-0.5, y1=len(df)+6.5,line=dict(color='Grey'))
  fig.add_shape(type="line",x0=-0.5, y0=len(df)+3.5, x1=len(mu_df)-0.5, y1=len(df)+3.5,line=dict(color='Grey'))      

  # Text Annotations
  fig.add_annotation(x=1, y=len(df)+8,text="B cell",xref='paper',yref='y',showarrow=False,xanchor='left')
  fig.add_annotation(x=1, y=len(df)+5,text="T cell (I)",xref='paper',yref='y',showarrow=False,xanchor='left')
  fig.add_annotation(x=1, y=len(df)+2,text="T cell (II)",xref='paper',yref='y',showarrow=False,xanchor='left')

  # Plot Epitope Regions
  epitope_df = stack_df[stack_df['epitope'] != 'None']
  for ind in epitope_df.index:
    if 'B cell' in stack_df['epitope'][ind]: 
      # B-Cell Epitope 
      fig.add_shape(type="rect",x0=stack_df['X_value'][ind]-0.5, 
              y0=len(df)+6.5, x1=stack_df['X_value'][ind]+0.5, y1=len(df)+9.5,
              line=dict(width=0),fillcolor="Grey")    
    if 'T cell HLA class I' in stack_df['epitope'][ind]:  
      # T-Cell Epitope (MHC Class I)
      fig.add_shape(type="rect",x0=stack_df['X_value'][ind]-0.5, 
              y0=len(df)+3.5, x1=stack_df['X_value'][ind]+0.5, y1=len(df)+6.5,
              line=dict(width=0),fillcolor="DarkGrey") 
    if 'T cell HLA class II' in stack_df['epitope'][ind]:
      # T-Cell Epitope (MHC Class II)
      fig.add_shape(type="rect",x0=stack_df['X_value'][ind]-0.5, 
              y0=len(df)+0.5, x1=stack_df['X_value'][ind]+0.5, y1=len(df)+3.5,
              line=dict(width=0),fillcolor="Black")       
       

  fig.update_layout(title='Mutation Profile per Sample')            
  fig.update_xaxes(title_text="Mutation")
  fig.update_yaxes(title_text="Virus Sample ID")          
  fig.write_html("output/07_mutation_profile.html")
  toc = time.perf_counter()
  print('Mutation Profile Viz Runtime: {:.4f} seconds'.format(toc-tic))



def epi_region(epitopes, epitope, protein, df):
  '''Returns 3 dataframes in the shape of df'''
  rdf = df.copy()
  id_df = df.copy()
  id_df[:] = ''
  epi_dict = {'B cell':0, 'T cell (class I)':1, 'T cell (class II)':2}
  for line in epitopes:
    if line[0] == epitope and line[1] == protein: 
      m = int(line[2])
      n = int(line[3])+1
      df.iloc[4:,m:n] = epitope + ' | '
      rdf.iloc[epi_dict[epitope],m:n] = line[4] # epitope restriction
      id_df.iloc[4:,m:n] = epitope + line[5] # B cell04
      
  epitope_region = df.mask(df!=epitope + ' | ', '')
  epitope_restriction = rdf.mask(df.astype('str').iloc[epi_dict[epitope]:epi_dict[epitope]+1]=='None','')
  epitope_id = id_df#.mask(df.astype('str')=='None','')
  
  return  epitope_region, epitope_restriction, epitope_id




def plot_mutations(fasta_DF, meta_DF, epitopes, proteins):
  ''' Creates a Plotly heatmap of the mutations '''

  
  display_size = 50
  
  for nprot in range(1,len(proteins) + 1): 
    prot_DF = fasta_DF[int(len(fasta_DF)/4) * (nprot-1) : int(len(fasta_DF)/4) * nprot]
    header = [None,None,None,None]
    prot_DF['Collection date'] = header.extend(meta_DF['Collection date'].tolist())
    prot_DF[4:] = prot_DF[4:].sort_values(by='Collection date', ascending=False)
    meta_DF = meta_DF.sort_values(by='Collection date', ascending=False)
    
    for batch in range(1,int(np.ceil(len(prot_DF)/display_size)) + 1):      
      fasta_df = prot_DF[0:4].append(prot_DF[display_size * (batch - 1) + 4 : display_size * batch + 4])
      meta_df = meta_DF[display_size * batch - display_size : display_size * batch]
      fasta_df.reset_index(drop=True, inplace=True)
      meta_df.reset_index(drop=True, inplace=True)
      print('\nCreating visualizations for {} protein batch #{}...'.format(fasta_df['name'][3], batch))

    
      df = pd.DataFrame(fasta_df['sequence'].str.replace('=','X'))
      #df = pd.DataFrame(df['sequence'].str.replace('-','X'))
      df = pd.DataFrame(df['sequence'].str.replace('_','X'))

    
      # Convert sequence strings into values of a list of lists
      aa = []
      for i in range(len(df)):
        aa.append(list(df['sequence'][i]))
    
  
      z = copy.deepcopy(aa)             
      aa = pd.DataFrame(aa)
      aa_code = aa.mask(aa=='W', '')
  
      # amino acids and other symbols to decimal number dictionary
      # equivalent to the rainbow colorscale
      num_dict = {  
                  'A':0.25, 'C':0.01, 'D':0.05, 'E':0.22, 'F':0.35, 
                  'G':0.80, 'H':0.52, 'I':0.50, 'K':0.95, 'L':0.40,
                  'M':0.18, 'N':0.60, 'P':0.75, 'Q':0.65, 'R':0.85,
                  'S':0.54, 'T':0.57, 'V':0.18, 'W':0.12, 'Y':0.20,
                  '-':0.30, 'X':None, 
                 }
  
      # Replace amino acids with integers
      for key in num_dict:
        replace_with(z,key,num_dict[key])
   
      # Create hover data 
      zdf = pd.DataFrame(z)            
      epitope_names = np.unique(np.array(epitopes)[:,0]).tolist()
      bcell_mask, bcell_res, bcell_id = epi_region(epitopes, epitope_names[0], fasta_df['name'][3], zdf)
      tcell1_mask, tcell1_res, tcell1_id = epi_region(epitopes, epitope_names[1], fasta_df['name'][3], zdf)
      tcell2_mask, tcell2_res, tcell2_id = epi_region(epitopes, epitope_names[2], fasta_df['name'][3], zdf)
      epitope_mask = bcell_mask + tcell1_mask + tcell2_mask
      # Epitope restriction mask
      restriction = bcell_res
      restriction.iloc[1] = tcell1_res.iloc[1]
      restriction.iloc[2] = tcell2_res.iloc[2]
       
      # mutation name hover data  e.g. N501Y  
      mutation = aa.iloc[3] + (aa.columns+1).astype(str) + aa    
      mutation.iloc[:4] = ''      
  
      # radical and conservative hover data
      radical_mask = aa.iloc[3] + aa # e.g. DA
      radical_mask = radical_mask.mask(radical_mask.isin(radical_list('grantham')), 'Radical')  
      radical_mask = radical_mask.mask(radical_mask!='Radical', 'Conservative')  
      radical_mask.iloc[:4] = ''
      
      # Color code the mutations Red (radical) and Green (conservative)
      radicality_color = zdf.mask(radical_mask=='Radical',0.95)
      radicality_color = radicality_color.mask(radical_mask=='Conservative',0.54)
      radicality_color = radicality_color.mask(aa=='X',None)
      zdf.iloc[4:] = radicality_color.iloc[4:] 
              
      # Create a Heatmap
      layout = go.Layout(
          yaxis=dict(dtick=1),
          xaxis=dict(range=[0, 150]))
  
      y = meta_df['Accession ID'].tolist()[::-1]
      y.extend([fasta_df['name'][3], 
                fasta_df['name'][2], 
                fasta_df['name'][1], 
                fasta_df['name'][0]]) 

  
      fig = go.Figure(data=go.Heatmap(
                x = np.arange(1,len(z[0])+1),
                y = y[::-1],
                z = zdf,
                customdata = np.dstack((aa_code, 
                                        bcell_mask, 
                                        tcell1_mask, 
                                        tcell2_mask, 
                                        restriction,
                                        mutation,
                                        radical_mask)), 
                type = 'heatmap',
                colorscale = 'rainbow',
                xgap = 0.5,
                ygap = 0.5,
                hoverongaps = False,
                showscale = False,
                hovertemplate = "%{y}<br>"
                                "%{x}<br>"
                                "%{customdata[0]}<br>"
                                "%{customdata[4]}<br>"
                                "%{customdata[1]}"
                                "%{customdata[2]}"
                                "%{customdata[3]}<br>"
                                "%{customdata[5]}<br>"
                                "%{customdata[6]}<br>"
                                "<extra></extra>"),
                layout = layout)
                 
      fig.update_layout(title='Sars-Cov-2 Mutations Viz | {} Protein | {} | {} to {}'.format(fasta_df['name'][3], 
                        meta_df['Location'][0].split('/')[1].rstrip().lstrip(), 
                        meta_df['Collection date'][len(meta_df)-1],
                        meta_df['Collection date'][4]),
                        xaxis=dict( rangeslider=dict(
                        visible=True)))    
                
      fig.update_xaxes(title_text="Site", rangeslider_thickness = 0.025)
      fig.update_yaxes(title_text="Accession ID") 
      fig['layout']['yaxis']['autorange'] = "reversed"
      fig.write_html("output/08_mutations/mutation_profile_{}_{}.html".format(fasta_df['name'][3],batch))         

      # Count of epitope radical hits
      bcell_rhits = np.zeros([len(bcell_id)], dtype='int')
      tcell1_rhits = np.zeros([len(tcell1_id)], dtype='int')
      tcell2_rhits = np.zeros([len(tcell2_id)], dtype='int')
      for row in range(4,len(bcell_id)):
        bcell_rhits[row] = len(bcell_id.mask(radical_mask!='Radical','').iloc[row].unique()) - 1
        tcell1_rhits[row] = len(tcell1_id.mask(radical_mask!='Radical','').iloc[row].unique()) - 1
        tcell2_rhits[row] = len(tcell2_id.mask(radical_mask!='Radical','').iloc[row].unique()) - 1
    
  
      # Create a summary table
      table_df = pd.DataFrame({ 'accession_id':y[::-1], 
                                'mutations':(aa!='X').sum(axis=1),
                                'mutations_epi':(aa.mask(epitope_mask=='','X') != 'X').sum(axis=1),
                                'r_mutations':(radical_mask=='Radical').sum(axis=1),
                                'bcell_r_mutations':bcell_rhits,
                                'tcell_r_mutations':tcell1_rhits + tcell2_rhits,
                                'epitope_r_mutations':bcell_rhits + tcell1_rhits + tcell2_rhits})
  
      summary_table(table_df, fasta_df, meta_df, batch)



  
def summary_table(table_df, fasta_df, meta_df, batch):
  
  header_values = [ 'Accession ID',
                    'Number of amino<br>acid substitutions<br>in the protein',
                    'Number of amino<br>acid substitutions<br>within epitope regions',
                    'Number of amino<br>acid substitutions<br>with radical hits',
                    'Number of B cell epitopes<br>with radical hits',
                    'Number of T cell epitopes<br>with radical hits',
                    'Total number of epitopes<br>with radical hits'
                  ] 
  
  table_plot = go.Figure(data=[go.Table(
      header=dict(values=header_values,
                  fill_color='royalblue',
                  align=['left','center'],
                  font=dict(color='white')),
      cells=dict(values=[
                        table_df[4:].accession_id, 
                        table_df[4:].mutations,
                        table_df[4:].mutations_epi, 
                        table_df[4:].r_mutations,
                        table_df[4:].bcell_r_mutations, 
                        table_df[4:].tcell_r_mutations, 
                        table_df[4:].epitope_r_mutations,
                        ],
                 fill_color=[['white','lightgrey','white','lightgrey']*len(table_df[4:])],
                 align=['left','center'])) ])  
  
  table_plot.update_layout(title='Summary Table | {} Protein | {} | {} to {}'.format(fasta_df['name'][3], 
                    meta_df['Location'][0].split('/')[1].rstrip().lstrip(),
                    meta_df['Collection date'][len(meta_df)-1],
                    meta_df['Collection date'][4]),
                    )
  
  table_plot.write_html('output/09_tables/table_{}_{}.html'.format(fasta_df['name'][3],batch))

  

def extract_region(meta_df):
  if len(meta_df.split('/')) > 2:
    region = meta_df.split('/')[2].lstrip().rstrip()
    
    # Enumerate all possible names of the regions
    ph_regions_dict = { 0  : ['NCR', 'Manila', 'National Capital Region'], 
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
    
    for key,value in ph_regions_dict.items():
      if region in value:
        return key     

  
  
def plot_variants(meta_DF):
  '''Creates a bar chart and a geographical map of the variants per region'''
  
  # Create a Regions column in the meta data dataframe
  meta_DF['Region_id'] = meta_DF['Location'].apply(extract_region)
  meta_df = meta_DF[~meta_DF['Region_id'].isna()]
  meta_df['Region_id'] = meta_df['Region_id'].astype(int)
  meta_df = meta_df.sort_values(by='Collection date', ascending=False)
  meta_df.reset_index(drop=True, inplace=True)
  meta_df['Collection date'] = meta_df['Collection date'].astype('datetime64')
  meta_df['Collection date'] = meta_df['Collection date'].dt.date

  # Variants of Concern and Variants of Interest 
  voc_voi = {   
                'B.1.1.7'  : 'VOC Alpha B.1.1.7 UK',
                'B.1.351'  : 'VOC Beta B.1.351 South Africa',
                'P.1'      : 'VOC Gamma P.1 Brazil',
                'B.1.617.2': 'VOC Delta B.1.617.2 India',
                'B.1.427'  : 'VOI Epsilon B.1.427 USA',
                'B.1.429'  : 'VOI Epsilon B.1.429 USA',
                'P.2'      : 'VOI Zeta P.2 Brazil',
                'B.1.525'  : 'VOI Eta B.1.526 UK/Nigeria',
                'P.3'      : 'VOI Theta P.3 Philippines',
                'B.1.526'  : 'VOI Iota B.1.526 USA',
                'B.1.617.1': 'VOI Kappa B.1.617.1 India',
                'C.37'     : 'VOI Lambda C.37 Peru',
                'Others'   : 'Others'
            }

  lineages_all = meta_df['Lineage'].unique()
  lineages = list(set(lineages_all).intersection(voc_voi.keys()))
  lineages.sort()
  meta_df = meta_df.replace(list(set(lineages_all).difference(lineages)), 'Others')

  ph_regions = [  'NCR', 'CAR', 'Region I', 'Region II', 'Region III', 
                  'Region IV-A', 'Region IV-B', 'Region V', 'Region VI', 'Region VII', 
                  'Region VIII', 'Region IX', 'Region X', 'Region XI', 'Region XII', 
                  'Caraga', 'BARMM']

  ph_regions_latitude = [ 14.5472655, 17.35971005, 15.855663, 16.7305959, 15.5, 
                          14.16052895, 13.072377, 13.42289175, 11, 10.6419531, 
                          11.9, 8.5, 8.3, 7.0428208, 6.2965755, 
                          9.2471392, 5.2]

  ph_regions_longitude = [121, 121.07525292, 120.1341247, 121.54067657, 121.13104248,
                         121.24256648, 121.3276166, 123.41103829, 122.58101501, 123.938142,
                         125,  123.4243754, 124.65677618, 125.580953, 124.9860759,
                         125.85578189, 120.02840346] 
         
  latitude = ph_regions_latitude
  longitude = ph_regions_longitude             
  regional_df = pd.DataFrame({'region':ph_regions,
                              'longitude':longitude,
                              'latitude':latitude})

  # Insert empty columns of unique lineages                            
  regional_df = regional_df.join(pd.DataFrame(0, index=range(len(regional_df)), columns=lineages_all))
  regional_df['Others'] = 0
  # Calculate the count per lineage per region
  groupby = pd.DataFrame(meta_df.groupby(['Lineage', 'Region_id'])['Accession ID'].count())
  # Insert the resulting counts into the regional_df
  for i in range(len(groupby)):
    lineage = groupby.index[i][0]
    region_id = groupby.index[i][1]
    regional_df[lineage][region_id] = groupby['Accession ID'][i]

  regional_df = regional_df[['region','longitude','latitude','Others']].join(regional_df[lineages])

  fig = make_subplots(rows=1, cols=2, 
                      specs=[[{'type':'scattergeo'}, {'type':'xy'}]])

  i = 0
  bar_data = [] 
  colors = ['red','blue','green','lightblue','magenta','yellow','lightgreen','gold','royalblue','maroon','whitesmoke'] 
  lineages.append('Others') 
   
  for lineage in lineages:    
    fig.add_trace(go.Scattergeo(
      lon = regional_df.longitude,
      lat = regional_df.latitude,
      text = regional_df.region, 
      name = voc_voi[lineage], 
      showlegend = False,
      customdata = regional_df[lineage],
      hovertemplate = 'Area: %{text}<br>'
                      'Occurence: %{customdata}',
      legendgroup = 'group{}'.format(i),
      marker = dict(
        size = regional_df[lineage]*50/322+2,
        color = colors[i], 
        line_width = 0,
        opacity  = 0.5,
      )),
      row=1, col=1)
      
    fig.add_trace(go.Bar( 
      name = voc_voi[lineage], 
      x = regional_df.region, 
      y = regional_df[lineage], 
      legendgroup = 'group{}'.format(i),
      marker_color = colors[i]),
      row=1, col=2)
      
    meta_df = meta_df.replace(lineage, voc_voi[lineage])
    i = i + 1    

  meta_df = meta_df.rename(columns={'Lineage':'Variant'})
  fig.update_layout(title = 'Sars-Cov-2 Variants per Region | {} | {} to {} | {} GISAID Sequences'.format(meta_df['Location'][0].split('/')[1].rstrip().lstrip(),
                          meta_df['Collection date'][len(meta_df)-1],
                          meta_df['Collection date'][0],
                          len(meta_DF)),
                    barmode = 'stack',
                    geo_resolution = 50,
                    geo_scope = 'asia',
                    geo_showframe = True,
                    geo_showcoastlines = True,
                    geo_landcolor = "rgb(217, 217, 217)",
                    geo_countrycolor = "black" ,
                    geo_coastlinecolor = "black",
                    geo_projection_type = 'mercator',
                    geo_lonaxis_range = [ 117.0, 127.0 ],
                    geo_lataxis_range = [ 5.0, 20.0 ],
                    legend = dict(
                        yanchor = "top",
                        y = 0.99,
                        xanchor = "left",
                        x = 0.8
                    ),
                    legend_title_text = 'Variants',
                    )
                            
  fig['layout']['yaxis']['title']='Count'
  fig.write_html('output/10_geoplot_variants.html') 
  
  variants_per_region(meta_df, ph_regions, colors)
  
  
  
def variants_per_region(meta_df, regions, colors):
  '''Create an area plot of the variants over time for each region'''
  
  plt.style.use('seaborn')
  groupby = meta_df[[ 'Variant',
                      'Collection date',
                      'Region_id',
                      'Accession ID']].groupby(by=[ 'Collection date',
                                                    'Variant',
                                                    'Region_id'], as_index=False).count()  
                                                 
  i = 0                                                  
  for region in regions:
    reg = groupby[groupby['Region_id']==i] 
    
    # If there's data
    if len(reg) > 0:
      reg_pivot = reg.pivot(index=reg['Collection date'], columns='Variant')['Accession ID']
      reg_pivot = reg_pivot.replace(np.nan, 0)
      total_df = reg_pivot.sum(axis=1)
      # Normalize the occurences
      for variant in reg_pivot.columns:
        reg_pivot[variant] = reg_pivot[variant]/total_df
        
        
      # Put the 'Others' column at the end
      cols = reg_pivot.columns.tolist()
      cols.sort(reverse=True)
      reg_pivot = reg_pivot[cols]
      
  
      # Smoothen curves
      df = pd.DataFrame()
      new_index = np.arange((reg_pivot.index[-1] - reg_pivot.index[0]).days + 1)
      date_index = pd.date_range(start = reg_pivot.index[0], 
                      end = reg_pivot.index[-1],freq='D').date

      for variant in reg_pivot.columns:
        func = interp1d((reg_pivot.index - reg_pivot.index[0]).days, reg_pivot[variant], kind='cubic')
        df[variant] = func(new_index)
      
      df[df < 0] = 0     
      df.index = date_index
      reg_pivot = df

      
    # If no data  
    else:  
      reg_pivot = pd.DataFrame({'Collection date':['2020-10','2021-01','2021-03'],'No data':[0,0,0]})
      reg_pivot['Collection date'] = reg_pivot['Collection date'].astype('datetime64')
      reg_pivot['Collection date'] = reg_pivot['Collection date'].dt.date

    reg_pivot.plot(kind = 'area', 
                  colors = colors[len(colors)-len(reg_pivot.columns):])
                    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), 
                                  key=lambda t: t[0], 
                                  reverse=True))
    ax.legend(handles, labels)
    plt.ylim(0,1)
    plt.xlim(meta_df['Collection date'].min(), meta_df['Collection date'].max())
    plt.xlabel('Collection date')
    plt.legend(loc = 'upper left') 
    plt.title('{}'.format(region), 
                    fontsize=17, 
                    fontstyle='oblique') 
    plt.savefig('output/11_regions/{}_{}.png'.format(i,region)) #dpi = 1200
    i = i + 1
    
  
  