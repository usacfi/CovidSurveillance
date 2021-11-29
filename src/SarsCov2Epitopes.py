from classes.Sequence import Sequence
import src.SarsCov2Variants as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import os
import time
import copy
import warnings

warnings.filterwarnings('ignore')


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










def calculate_protein_diffs(seq_dict, gene_name, ref_gene, start_loc, df):
    for name in seq_dict:
        # trigger protein diff calculation if sequence is not reference
        seq = seq_dict[name]
        df = seq.calculate_protein_diffs(name, gene_name, ref_gene, name==gene_name, start_loc, df)
    return df












def combine_fasta_files(directory):
    
    os.system(f'rm {directory}/combined.fasta')
    
    count = 0
    for dirname, _, filenames in os.walk(directory):
        for filename in filenames:
            if 'sequences.fasta' in filename:
                # copy to outer folder
                os.system(f'cp {dirname}/{filename} {directory}/temp_fasta{count}.fasta')                               
                count = count + 1

    if count == 0:
        print('\nError: No data was found. Check the spelling of filenames and check if the files are not empty.\n')
        exit()
    else:    
        os.system(f'awk 1 {directory}/*.fasta > {directory}/combined.fasta')
        os.system(f'find {directory}/ -type f -name "temp_fasta*" -delete')
        print(f'\nFound {count} "sequences.fasta" files and combined them into 1 file "combined.fasta"\n')





    





def create_df(df, genes, epitopes, division_column='agg_division'):
    aa_diffs = []        
    for site in df.columns.values:
        if isinstance(site, int):
            for mutation in df[site].dropna().unique():
                radical_replacement = is_radical(mutation)
                if division_column in df.columns: region_instances = instances_per_region(df, site, mutation) 
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
                if division_column in df.columns: cols.extend(region_instances)
                aa_diffs.append(cols)
                
    ph_regions = ['NCR', 'CAR', 'Region I', 'Region II', 'Region III', 'Region IV-A',
                'Region IV-B', 'Region V', 'Region VI', 'Region VII', 'Region VIII',
                'Region IX', 'Region X', 'Region XI', 'Region XII', 'Caraga', 'BARMM']              
    column_names = ['site', 'mutation','instances', 'gene', 'epitope', 'replacement']
  
    if division_column in df.columns: column_names.extend(ph_regions)              
    new_df = pd.DataFrame(aa_diffs, columns=column_names)
    new_df = new_df.sort_values(by='site')
    new_df = new_df.reset_index(drop=True)

    return new_df











def create_sequences(seq_dict, start_loc, stop_loc):
    for name in seq_dict:
        seq = seq_dict[name]
        seq_dict[name] = Sequence(name, seq, start_loc, stop_loc) # instantiate Sequence












def fasta_to_heatmap(fasta_file, genes=['Spike','Envelope','Membrane','Nucleocapsid'], gene=0, 
                    meta_df=scv.combine_metadata(f'references/Sequences/Philippines/Epitope_Surveillance')):
                    
    '''
    Description: Converts a fasta_df into a dataframe with 3 columns as inputs to a heatmap
    
    
    Returns a heatmap dataframe with the following columns
    ID: virus isolate identifying string
    Site: location of the amino acid in the gene
    Amino Acid: 
    Radical: type of mutation whether conservative or radical
    Color: heatmap color
    Mutation: amino acid replacements with the format [original AA][site][new AA] e.g. 'N501Y'
    Variant: WHO labels or variants of concern and variants of interest
                    
    '''
    
 
    
    num_genes = len(genes)
    fasta_df = fasta_to_df(fasta_file)
    len_fasta = len(fasta_df)//num_genes
    meta_df = scv.lineage_to_variant(meta_df)
    #meta_df = scv.aggregate_divisions(meta_df)   
    meta_df['date'] = pd.to_datetime(meta_df['date']).dt.date
    meta_df.sort_values(by='date', ascending=False, inplace=True)
    meta_df.reset_index(drop=True, inplace=True)
    meta_df['batch'] = 'Header'    
    
    # Add a Batch column to limit the number rows shown 
    for batch in range(len_fasta//30 + 1):
        start = 30*batch
        end = 30*(batch+1)
        end_date = meta_df.date.iloc[start]   
        if end > len_fasta:
            start_date = meta_df.date.iloc[-1]
            print(start_date)
        else: start_date = meta_df.date.iloc[end]   
        meta_df.loc[start : end, ['batch']] = f'Batch {batch} ({start_date} - {end_date})' 
        
    meta_df = meta_df[['strain','variant','batch']]    
    fasta_df = fasta_df[int(gene * len_fasta) : int((gene + 1) * len_fasta)]  
    len_seq = len(fasta_df['sequence'].iloc[0])
    fasta_df = fasta_df.merge(meta_df, left_on='name', right_on='strain', how='left')    
    fasta_df.loc[0:3, ['batch']] = 'Header'


    
    # Create the initial lists
    rows = []
    columns = []
    values = []
    references = []
    variants = []
    batches = []
    for idx in range(len(fasta_df)):
        rows.extend([fasta_df['name'].iloc[idx]] * len_seq)
        columns.extend(list(range(1,len_seq+1)))
        values.extend(list(fasta_df['sequence'].iloc[idx]))
        references.extend(list(fasta_df['sequence'].iloc[3]))
        variants.extend([fasta_df['variant'].astype(str).iloc[idx]] * len_seq)
        batches.extend([fasta_df['batch'].astype(str).iloc[idx]] * len_seq)
        
    # Create the initial heatmap_df  
    heatmap_df = pd.DataFrame({
                                'ID': rows,
                                'Site': columns,
                                'Amino Acid': values,
                                'Reference': references,
                                'Variant': variants,
                                'Batch': batches,
                            })
                            
    heatmap_df['Amino Acid'] = heatmap_df['Amino Acid'].astype(str).str.replace('=','')
    heatmap_df['Amino Acid'] = heatmap_df['Amino Acid'].astype(str).str.replace('_','')
                 
    
    # Add a radical column
    heatmap_df['Replacement'] = heatmap_df['Reference'] + heatmap_df['Amino Acid']
    grantham_criteria = radical_list('grantham')
    heatmap_df['Radical'] = list(map(lambda x: x in grantham_criteria, heatmap_df['Replacement'].tolist()))
    heatmap_df['Radical'] = ['Radical' if x==True else 'Conservative' for x in heatmap_df.Radical.tolist()]
          
    # Add a color column
    heatmap_df['Color'] = heatmap_df['Radical'].tolist()
    
    # Create HLA restriction column
    
    
    # Add mutation column
    heatmap_df['Mutation'] = heatmap_df['Reference'] +  heatmap_df['Site'].astype(str) + heatmap_df['Amino Acid']
    
    # Further cleanup for the reference genes and the epitopes
    for gene in genes:
        heatmap_df.loc[heatmap_df.ID==gene, ['Radical']] = ''
        heatmap_df.loc[heatmap_df.ID==gene, ['Color']] = heatmap_df[heatmap_df.ID==gene]['Amino Acid'].tolist()
        heatmap_df.loc[heatmap_df.ID==gene, ['Mutation']] = ''
        heatmap_df.loc[heatmap_df.ID==gene, ['Variant']] = ''
        
    for epitope in ['B-Cell', 'T-Cell (class I)', 'T-Cell (class II)']:
        heatmap_df.loc[(heatmap_df.ID==epitope) & (heatmap_df['Amino Acid']=='W'), ['Amino Acid']] = '*'
        heatmap_df.loc[heatmap_df.ID==epitope, ['Radical']] = ''
        heatmap_df.loc[heatmap_df.ID==epitope, ['Color']] = '*'
        heatmap_df.loc[heatmap_df.ID==epitope, ['Mutation']] = ''
        heatmap_df.loc[heatmap_df.ID==epitope, ['Variant']] = ''
    
    return heatmap_df[['ID', 'Site', 'Amino Acid', 'Radical', 'Color', 'Mutation', 'Variant', 'Batch']]
    











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











def epi_region(epitopes, epitope, protein, df):
    '''Returns 3 dataframes in the shape of df'''
    rdf = df.copy()
    id_df = df.copy()
    id_df[:] = ''
    epi_dict = {'B cell':0, 
                'T cell (class I)':1, 
                'T cell (class II)':2
                }
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
  
  
  
  
  
  
  
  
  
  
  
  
def instances_per_region(df, site, mutation, division_column='agg_division'):
  
    instance_list = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]              
    location = df[df[site]==mutation][division_column]

    # Enumerate all possible names of the regions
    region_dict ={0  : ['NCR', 'Manila', 'National Capital Region'], 
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
  
  
  
  
  
  
  
  
  
    
    
    
def is_radical(mutation):
    ''' Classifies a mutation into 1 (radical) or 0 (conservative) '''

    replacement = mutation[0] + mutation[-1]
    reverse = mutation[-1] + mutation[0]

    # Based on Grantham's mean chemical difference index (1974)
    grantham=['RS', 'LS', 'LR', 'PR', 'AR', 'VR', 'GR', 'GL', 'GV', 'IS', 'IG', 'FS',
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










# Not used
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
  
  











def plot_mutations(fasta_DF, meta_DF, epitopes, proteins, date_column='date', id_column='gisaid_epi_isl',
                    division_column='agg_division'):
                    
  ''' Creates a Plotly heatmap of the mutations '''
  
  display_size = 50
  
  for nprot in range(1,len(proteins) + 1): 
    prot_DF = fasta_DF[int(len(fasta_DF)/4) * (nprot-1) : int(len(fasta_DF)/4) * nprot]
    header = [None,None,None,None]
    prot_DF[date_column] = header.extend(meta_DF[date_column].tolist())
    prot_DF[4:] = prot_DF[4:].sort_values(by=date_column, ascending=False)
    meta_DF = meta_DF.sort_values(by=date_column, ascending=False)
    
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
  
      y = meta_df[id_column].tolist()[::-1]
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
                        meta_df[division_column][0], 
                        meta_df[date_column][len(meta_df)-1],
                        meta_df[date_column][4]),
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












def radical_list(criteria='grantham'):
  
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
      
    
     
     
     
     
     
     
     
     
     
     
     
def replace_at(string, index, rep):   
  ''' Replace a character on a specific index of a string ''' 
  return string[0:index] + rep + string[index+1:]  
    














def replace_with(my_lists, old, new):
  ''' Replace a character in a list of lists or numpy array'''
  for i in range(len(my_lists)):
    my_lists[i] = [new if x==old else x for x in my_lists[i]]
    
  return my_lists
    
    
    
    
    
    
    
    
    
    
    
    
# Not used 
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


  
  
  
  
  
  
  
  
  
  
  

def summary_table(table_df, fasta_df, meta_df, batch, date_column='date', division_column='agg_division'):
  
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
        cells=dict( values=[
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
                    meta_df[division_column][0].split('/')[1].rstrip().lstrip(),
                    meta_df[date_column][len(meta_df)-1],
                    meta_df[date_column][4]),
                    )

    table_plot.write_html('output/09_tables/table_{}_{}.html'.format(fasta_df['name'][3],batch))





  
  
  
  
  
  
  
  
  
def which_protein(site, proteins):
    ''' proteins is a list of lists, e.g.
      [['Surface',123,456],
       ['Envelope',432,765],
       ['Membrane',234,567]] '''
  
    for protein in proteins:
        if site in range(int(protein[1]), int(protein[2])):
            return protein[0]   
    raise ValueError(f'Site location {site} is outside the (nucleotide) regions of the relevant proteins')     
  
  
  
  




  


