
# 1/5/21 version
# parse_alignments_v2.py
# parse info about alignments
# works with Python 3.7

# ----------------------------------------------------------------------
import argparse

parser = argparse.ArgumentParser(description="parses SNPs from alignments")

parser.add_argument('-aln', '--alignments_fasta', dest='aln', help='input alignment file', required=True)

parser.add_argument('-o', '--output', dest='out', help='base name for output files', required=True)

parser.add_argument('-ref', '--reference_name', dest='ref',
                    help='name of reference sequence in alignment to compare with other sequences; '
                         'otherwise uses first sequence in alignment as reference')

parser.add_argument('-eps', '--epitope_regions', dest='eps', help='text file that lists epitopes to check')

args = parser.parse_args()

# ----------------------------------------------------------------------

#######################################################################

class Sequence:
    ref_name = None
    rprot = None

    def __init__(self, name, sequence):
        # print("constructing amplicon: " + name)
        self.name = name
        self.sequence = sequence.upper()
        self.seq_no_gaps = remove_gaps(self.sequence)
        # TODO: check that insertions and deletions are handled correctly
        self.protein = translate(self.sequence)
        self.prot_diffs = []
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

        output = [self.name, self.protein, ','.join(self.prot_diffs), ','.join(e_positions), ','.join(e_seqs),
                  ','.join(e_diffs)]
        return '\t'.join(map(str, output))

    def add_epitope(self, start, end):
        new_epi = Epitope(self, start, end)
        self.epitopes.append(new_epi)

class Epitope:
    def __init__(self, seq, start, end):
        self.start = int(start)
        self.end = int(end)
        self.sequence = None
        self.diffs = []
        self.find_changes(self.start, self.end, seq.protein)

    def find_changes(self, start, end, prot):
        rprot = Sequence.rprot
        count = 0
        seq = ''
        for i in range(0, end):
            if count < (start-1):
                if rprot[i] != '-':
                    count += 1
                continue
            seq += prot[i]
            if rprot[i] != '-':
                count += 1
            if rprot[i] != prot[i]:
                diff = rprot[i] + str(count) + prot[i]
                # print(diff)
                self.diffs.append(diff)
        self.sequence = seq

def remove_gaps(seq):
    no_gaps = seq.replace("-", "")
    return no_gaps

def translate(seq):

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
        #TODO: add wobble codons
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table:
                protein += table[codon]
            else:
                protein += 'X'
    else:
        print("ERROR: seq incorrect length for translating" + " --length: " + str(len(seq)))

    return protein

def parse_alignment(fasta_file):
    seq_dict = dict()
    with open(fasta_file) as data:
        seq = ''
        name = ''
        for line in data:
            line = line.rstrip()
            if line.startswith('>'):
                if Sequence.ref_name is None:
                    Sequence.ref_name = line[1:] # first sequence is the reference

                if seq != '':
                    new_seq = Sequence(name, seq)
                    seq_dict[name] = new_seq
                # start creating next Sequence object:
                name = line[1:]
                seq = ''

            else:
                line = ''.join(line.split()) # remove whitespace
                seq += line.upper() # add to seq, which may be on multiple lines

        # add last seq:
        new_seq = Sequence(name, seq)
        seq_dict[name] = new_seq

    return seq_dict

def parse_protein_diffs(seq_dict):
    rprot = Sequence.rprot
    for key in seq_dict:
        seq = seq_dict[key]

        if seq is not ref_seq:
            prot = seq.protein
            count = 0
            for i in range(0,len(rprot)-1):
                if rprot[i] != '-':
                    count += 1
                if rprot[i] != prot[i]:
                    diff = rprot[i] + str(count) + prot[i]
                    seq.prot_diffs.append(diff)

def parse_epitopes(epitope_file_name, seq_dict):
    with open(epitope_file_name) as eps:
        for line in eps:
            line = line.rstrip()
            values = line.split('-')
            start = values[0]
            end = values[1]

            for seq in seq_dict.values():
                seq.add_epitope(start, end)

#######################################################################

if __name__ == "__main__":

    seq_dict = parse_alignment(args.aln)
    print("reference name: " + Sequence.ref_name)
    ref_seq = seq_dict[Sequence.ref_name]
    Sequence.rprot = ref_seq.protein

    if args.eps is not None:
        parse_epitopes(args.eps, seq_dict)

    parse_protein_diffs(seq_dict)

    with open(args.out + '.xls', 'w') as out:
        for seq in seq_dict.values():
            out.write(str(seq) + '\n')
            #print(str(seq))

    #TODO: output separate file with epitope information; different format?
