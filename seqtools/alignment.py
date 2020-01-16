from Bio import AlignIO
import os
from shutil import copyfile
import subprocess
import tempfile
import uuid


TMP_DIR = tempfile.mkdtemp(dir='/lrlhps/scratch/c259188/tmp/')


def run_clustalw2(input_seqs, output_file=None, tmp_dir=TMP_DIR):
    """
    Multiple Sequence Alignment
    Run ClustalW2 on commandline
    input_seqs is a list of lists; each item is [name, sequence]
    Alignment object is returned
    output_file is written if provided
    """
    filename = uuid.uuid4().hex[:10]
    fa_file = os.path.join(tmp_dir, '{0}.fa'.format(filename))
    with open(fa_file, 'w') as input_fh:
        for seq in input_seqs:
            input_fh.write('>{0}\n{1}\n'.format(
                seq[0], seq[1],
            ))
    subprocess.check_call([
        'clustalw2', '-infile={0}'.format(fa_file),
        '-outorder=input', 
    ])
    os.remove(fa_file)
    aln_file = os.path.join(tmp_dir, '{0}.aln'.format(filename))
    alignment = AlignIO.read(aln_file, 'clustal')
    if output_file:
        copyfile(aln_file, output_file)
    os.remove(aln_file)
    return alignment


def run_needle(input_seqs,  output_file=None, tmp_dir=TMP_DIR):
    """
    Pairwise Sequence Alignment
    Run Emboss Needle on commandline
    input_seqs is a list of 2 sequences; each item is [name, sequence]
    Alignment object is returned
    output_file is written if provided
    """
    if not len(input_seqs) == 2:
        raise IOError('2 input sequences required, ' \
                      'if more consider running clustalw2 ' \
                      'multiple sequence alignment')
    filename = uuid.uuid4().hex[:10]
    fa_file1 = os.path.join(tmp_dir, '{0}.fa'.format(filename))
    filename = uuid.uuid4().hex[:10]
    fa_file2 = os.path.join(tmp_dir, '{0}.fa'.format(filename))
    filename = uuid.uuid4().hex[:10]
    outfile = os.path.join(tmp_dir, '{0}.needle'.format(filename))
    with open(fa_file1, 'w') as input_fh:
        input_fh.write('>{0}\n{1}\n'.format(
            input_seqs[0][0], input_seqs[0][1],
        ))
    with open(fa_file2, 'w') as input_fh:
        input_fh.write('>{0}\n{1}\n'.format(
            input_seqs[1][0], input_seqs[1][1],
        ))
    subprocess.check_call([
        'needle', '-asequence={0}'.format(fa_file1),
        '-bsequence={0}'.format(fa_file2),
        '-gapopen=10', '-gapextend=0.5',
        '-outfile={0}'.format(outfile)
    ])
    os.remove(fa_file1); os.remove(fa_file2)
    alignment = AlignIO.read(outfile, 'emboss')
    if output_file:
        copyfile(outfile, output_file)
    os.remove(outfile)
    return alignment


def get_contigs_from_clustal(alignment, csv_out, base_sequence=None):
    """
    alignment = Multiple Sequence Alignment Object 
      read in using AlignIO.read
    """
    contigs = []
    contig = ''
    match = False
    start = 0
    gaps = 0
    if base_sequence:
        base_seq = [i for i in alignment if i.id == base_sequence][0]
    else:
        base_seq = alignment[0]
    for num,nt in enumerate(base_seq):
        if nt == '-':
            gaps += 1
        if all([nt == i[num] for i in alignment]):
            if match:
                contig += nt
            else:
                contig = nt
                start = num + 1
                match = True
        else:
            if match:
                contigs.append({
                    'start_wrt_alignment': start-1,
                    'end_wrt_alignment': num,
                    'sequence': contig,
                    'length': len(contig),
                    'start_wrt_' + base_seq.id: start - gaps,
                    'end_wrt_' + base_seq.id: num - gaps,
                })
                match = False
    headers = ['start_wrt_alignment', 'end_wrt_alignment', 
               'start_wrt_' + base_seq.id, 'end_wrt_' + base_seq.id,
               'length', 'sequence']
    with open(csv_out, 'w') as output_fh:
        output_fh.write('{0}\n'.format(','.join(headers)))
        for contig in contigs:
            output_fh.write('{0}\n'.format(
                ','.join([str(contig[i]) for i in headers])
            ))


if __name__ == '__main__':
    pass
