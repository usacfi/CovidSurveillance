class Epitope:
  def __init__(self, sequence, start, end, ref_protein):
    self.start = start
    self.end = end
    self.sequence = None
    self.diffs = []

    protein = sequence.protein
    self.calculate_epitope_diffs(start, end, protein, ref_protein)

  def calculate_epitope_diffs(self, start, end, protein, ref_protein):
    count = 0
    seq = ''

    for i in range(0, end):
      if count < (start-1):
        if ref_protein[i] != '-':
          count += 1
        continue
      
      seq += protein[i]

      if ref_protein[i] != '-':
        count += 1

      if ref_protein[i] != protein[i]:
        diff = ref_protein[i] + str(count) + protein[i]
        self.diffs.append(diff)

    self.sequence = seq