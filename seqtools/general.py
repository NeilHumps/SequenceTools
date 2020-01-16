"""
Useful python scripts
"""


import sys


CODON_TO_AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*', 
}

SINGLE_TO_TRIPLET = {
    'G': 'Gly',
    'P': 'Pro',
    'A': 'Ala',
    'V': 'Val',
    'L': 'Leu',
    'I': 'Ile',
    'M': 'Met',
    'C': 'Cys',
    'F': 'Phe',
    'Y': 'Tyr',
    'W': 'Trp',
    'H': 'His',
    'K': 'Lys',
    'R': 'Arg',
    'Q': 'Gln',
    'N': 'Asn',
    'E': 'Glu',
    'D': 'Asp',
    'S': 'Ser',
    'T': 'Thr',
    '*': '*',
}


DNA_COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
    'M': 'K',
    'R': 'Y',
    'W': 'W',
    'S': 'S',
    'Y': 'R',
    'K': 'M',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'B': 'V',
    'X': 'X',
}

DNA_RNA = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'U',
}

HUMAN_CHROMOSOMES = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
    'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
    'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
    'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM',
]

HUMAN_CHROMOSOMES_ENS = [
    '1', '2', '3', '4', '5', '6',
    '7', '8', '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18',
    '19', '20', '21', '22', 'X', 'Y', 'MT',
]

MOUSE_CHROMOSOMES_ENS = [
    '1', '2', '3', '4', '5', '6',
    '7', '8', '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18',
    '19', 'X', 'Y', 'MT',
]


def load_csv(input_csv, sep=',', headers=None):
    """
    Load delimited format files into a list of dicts
     - input_csv = input filename
     - sep = column separator (default=',')
     - if not headers, first line will be used as header row,
       otherwise headers can be specified and must correspond to column contents
    """
    output_csv = []
    with open(input_csv, 'rU') as input_fh:
        for line in input_fh:
            if not headers:
                headers = line.strip('\n\r').split(sep)
            else:
                output_csv.append(dict(zip(headers, line.strip('\n\r').split(sep))))
    print(' - {0} lines parsed with {1} columns'.format(len(output_csv), len(headers)))
    return headers, output_csv


def write_csv(filename, dict_list, sep=',', headers=None, write_headers=True):
    """
    Write a list of dicts into delimited format
     - filename = output filename
     - dict_list = python list of dicts (headers must match key names)
     - sep = column separator (default=',')
     - if not headers, keys of first dict will be used as header row,
       otherwise headers can be specified and must exist in every dict
    """
    if not headers:
        headers = dict_list[0].keys()
    rowcount = 0
    colcount = len(headers)
    with open(filename, 'w') as output_fh:
        if write_headers:
            output_fh.write(sep.join(headers) + '\n')
        for line_dict in dict_list:
            output_fh.write(sep.join(map(str, [line_dict[i] for i in headers])) + '\n')
            rowcount += 1
    print(' - {0} lines x {1} columns written to:\n   - {2}'.format(
        rowcount, colcount, filename
    ))


def translate(input_sequence, out_format='single'):
    """
    Translate DNA sequence into single letter amino acid code
     - out_format can be single or triplet (default=single)
    """
    output_seq = []
    for pos in range(0, len(input_sequence), 3):
        codon = input_sequence.upper()[pos:pos+3]
        try:
            output_seq.append(CODON_TO_AA[codon])
        except KeyError:
            if len(codon) == 3:
                raise IOError('Invalid character in input sequence')
    if out_format == 'single':
        return ''.join(output_seq)
    elif out_format == 'triplet':
        return '-'.join([SINGLE_TO_TRIPLET[i] for i in output_seq])
    else:
        raise IOError('out_format must be single or triplet')


def transcribe(input_sequence):
    """
    Transcribe DNA sequence into RNA
    """
    output_seq = []
    for nt in input_sequence.upper():
        try:
            output_seq.append(DNA_RNA[nt])
        except KeyError:
            raise IOError('Untranscribable character in input sequence')
    return ''.join(output_seq)


def rc(input_sequence, complement=True, reverse=True):
    """
    Convert DNA sequence to reverse complement
    """
    output_seq = []
    if reverse:
        sequence = list(reversed([i for i in input_sequence.upper()]))
    else:
        sequence = [i for i in input_sequence.upper()]
    if complement:
        for nt in sequence:
            try:
                output_seq.append(DNA_COMPLEMENT[nt])
            except KeyError:
                raise IOError('Invalid DNA Character')
    else:
        output_seq = sequence
    return ''.join(output_seq)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'rc':
            print(rc(sys.argv[2]))
        elif sys.argv[1].lower() == 'translate':
            print(translate(sys.argv[2]))
        elif sys.argv[1].lower() == 'transcribe':
            print(transcribe(sys.argv[2]))
        else:
            print('Invalid input function, rc, translate or transcribe are valid')
    else:
        print('Invalid input function, rc, translate or transcribe are valid')
