from classes.Sequence import Sequence
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import folium
from folium.plugins import MarkerCluster
import osmnx as ox


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



def create_df(df, genes, epitopes):
  aa_diffs = []        
  for site in df.columns.values:
    if isinstance(site, int):
      for mutation in df[site].dropna().unique():
        radical_replacement = is_radical(mutation)
        region_instances = instances_per_region(df, site, mutation) 
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
        cols.extend(region_instances)
        aa_diffs.append(cols)
                
  new_df = pd.DataFrame(aa_diffs, columns=['site', 'mutation','instances', 'gene', 'epitope', 'replacement',
                                            'NCR', 'CAR', 'Region I', 'Region II', 'Region III', 'Region IV-A',
                                            'Region IV-B', 'Region V', 'Region VI', 'Region VII', 'Region VIII',
                                            'Region IX', 'Region X', 'Region XI', 'Region XII', 'Caraga', 'BARMM'])
  new_df = new_df.sort_values(by='site')
  new_df = new_df.reset_index(drop=True)
  return new_df
  
  
  
def geoplot_mutations(df):
  # https://www.openstreetmap.org/export
  llat = 5.091
  ulat = 21.412
  llon = 116.499
  ulon = 127.134
  center = [(llat+ulat)/2,(llon+ulon)/2]
  
  regions = df.columns[6:].values
  radical = np.zeros([len(regions)], dtype='int')
  conservative = np.zeros([len(regions)], dtype='int')
  vaccine_index = np.zeros(len(regions))
  #latitude = np.zeros([len(regions)])
  #longitude = np.zeros([len(regions)])
  
  latitude = [14.5472655, 17.35971005, 15.855663, 16.7305959, 15.1954795, 14.16052895,
             13.072377, 13.42289175, 10.69777705, 10.6419531, 11.8029647, 8.6549449,
             8.486114, 7.0428208, 6.2965755, 9.2471392, 5.8506231]
  
  longitude = [120.98448771, 121.07525292, 120.1341247, 121.54067657, 121.13104248,
               121.24256648, 121.3276166, 123.41103829, 122.58101501, 123.938142,
               124.8663766,  123.4243754, 124.65677618, 125.580953, 124.9860759,
               125.85578189, 120.02840346] 
  
  for i in range(len(regions)):
    radical[i] = df[df['replacement']=='Radical'][regions[i]].agg(sum)
    conservative[i] = df[df['replacement']=='Conservative'][regions[i]].agg(sum)
    vaccine_index[i] = radical[i]/(radical[i] + conservative[i] + 1e-9) 
    #latitude[i],longitude[i] = ox.geocode('Philippines, '+ str(regions[i]))
    
  regional_df = pd.DataFrame({'region':regions,
                              'longitude':longitude,
                              'latitude':latitude,
                              'radical':radical,
                              'conservative':conservative,
                              'index':vaccine_index})
  
  my_map = folium.Map(location = center, 
                      zoom_start=6, 
                      tiles='Stamen Terrain')
                      
  for j in range(len(regions)):
    if int(regional_df['radical'][j]) != 0:
      folium.CircleMarker([ regional_df['latitude'][j], regional_df['longitude'][j]], 
                        radius=regional_df['index'][j] * 50, 
                        color='red', 
                        opacity=0.5,
                        tooltip=regional_df['region'][j], 
                        popup='Vaccine Index: {:.4}'.format(regional_df['index'][j]),
                        fill=True, 
                        fillColor='red').add_to(my_map)
  
  my_map.save('output/08_mutation_geomap.html') 