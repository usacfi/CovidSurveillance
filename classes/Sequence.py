from classes.Epitope import Epitope

class Sequence:
  def __init__(self, name, sequence):
    self.name = name
    self.sequence = sequence
    self.sequence_no_gaps = self.remove_gaps(sequence)
    self.protein = self.translate(sequence)
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
      self.protein, 
      ','.join(self.protein_diffs), 
      ','.join(e_positions), 
      ','.join(e_seqs),
      ','.join(e_diffs)
    ]

    return '\t'.join(map(str, output))

  def remove_gaps(self, sequence):
    return sequence.replace("-", "")

  def translate(self, sequence):
    # TODO: add wobble codons
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
    # intantiate a new Epitope, also calculated epitope protein difference
    new_epi = Epitope(self, start, end, ref_protein)
    self.epitopes.append(new_epi)

  def calculate_protein_diffs(self, ref_name, ref_protein):
    count = 0

    for i in range(0, len(ref_protein)-1):
      if ref_protein[i] != '-':
        count += 1
      
      if ref_protein[i] != self.protein[i]:
        diff = ref_protein[i] + str(count) + self.protein[i]
        self.protein_diffs.append(diff)