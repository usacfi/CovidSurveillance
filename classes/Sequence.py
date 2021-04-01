from classes.Epitope import Epitope
from itertools import chain

class Sequence:
  def __init__(self, name, sequence, start_loc, stop_loc):
    self.name = name
    self.sequence = sequence
    # Trim sequence
    self.sequence = sequence[start_loc:stop_loc]
    with open('output/02_trimmed.fasta','a') as out:
      out.write('>'+ str(name) + '\n')
      out.write(str(self.sequence).rjust(stop_loc, '_') + '\n')
    self.sequence_no_gaps = self.remove_gaps(self.sequence) # not needed
    self.protein = self.translate(self.sequence)
    with open('output/03_protein.fasta','a') as out:
      out.write('>'+ str(name) + '\n')
      out.write(str(self.protein) + '\n')
    self.mutations = []
    self.protein_diffs = []
    self.epitopes = []

  def __str__(self):
    e_positions = []
    e_seqs = []
    e_diffs = []

    for e in self.epitopes:
      position = str(e.start) + '-' + str(e.end)
      e_positions.append(position)
      e_seqs.append(e.sequence)
      e_diffs.append(';'.join(e.diffs))

    output = [
      self.name, 
      #self.protein, # As discussed, this this not important as output
      ''.join(self.mutations), # technically exactly the same with protein_diffs, but with format
      #','.join(self.protein_diffs), # As discussed, this this not important as output
      ','.join(e_positions), 
      ','.join(e_seqs),
      ','.join(e_diffs)
    ]

    return '\t'.join(map(str, output))

  def remove_gaps(self, sequence):
    return sequence.replace("-", "")

  def translate(self, sequence):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        '---': '-',
        # wobble codons - based on codon chart and amino acid codes
        'TCN': 'S', 'CTN': 'L', 'CCN': 'P', 'CGN': 'R',
        'ACN': 'T', 'GTN': 'V', 'GCN': 'A', 'GGN': 'G',
    }

    protein = ''

    if len(sequence) % 3 == 0:
      # iterate through each codon
      for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += table[codon] if codon in table else 'X'
    else:
      print("ERROR: seq incorrect length for translating" + " --length: " + str(len(sequence)))
    
    return protein

  def add_epitope(self, start, end, ref_protein):
    # instantiate a new Epitope, also calculated epitope protein difference
    new_epi = Epitope(self, start, end, ref_protein)
    self.epitopes.append(new_epi)

  def calculate_protein_diffs(self, name, ref_name, ref_protein, is_ref, start_loc, df):
    count = 0

    for i in range(0, len(ref_protein)-1):
      if ref_protein[i] != '-':
        count += 1
      
      # if no difference and not the reference sequence
      if ref_protein[i] == self.protein[i] and not(is_ref) :
        self.mutations.append('-')

      # if proteins are different or if it is the reference sequence
      # technically there should be no difference, but based on the requested output,
      # we would like to see the reference sequence on top of the differences
      else:
        self.mutations.append(self.protein[i])
        # Assuming that 'X, -, _' are not a real mutations
        if self.protein[i] not in ['X', '-', '_'] and not(is_ref):
          protein_diff = ref_protein[i] + str(i+1) + self.protein[i]
          df.loc[int(df[df['ID']==name].index[0]), start_loc + i + 1] = protein_diff
          #print('{} | {} | {}'.format(int(df[df['ID']==name].index[0]), start_loc + i + 1, protein_diff))
          
    return df
        
    