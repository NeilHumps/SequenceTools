from Bio import SeqIO
import gzip
import os
import sys

from seqtools.general import rc, translate


def parse_gtf_dict(gtf_str):
    return {i.split(' "')[0]:i.split(' "')[1] for i in gtf_str.split('"; ')}


def gtf_to_gene_length(gtf_file, outfile, sum_type='transcript_longest'):
    """
    Get Gene length, 3 options:
    1. transcript_longest: Length of longest transcript (default)
    2. gene_unique: Length of all unique 'expressed DNA' i.e. count every base that is expressed once
    3. gene_all: Total length of all 'expressed DNA' i.e. if transcripts contain overlapping exons then count exon length multiple times
    """
    if not sum_type == 'transcript_longest':
        gene_dict = {}
        gene_count = 0
        with open(gtf_file) as input_fh:
            for line in input_fh:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'exon':
                        gtf_dict = parse_gtf_dict(fields[8])
                        if not gtf_dict['gene_name'] in gene_dict:
                            gene_count += 1
                            if gene_count % 1000 == 0:
                                print(' - {0} genes processed'.format(gene_count))
                            gene_dict[gtf_dict['gene_name']] = []
                        gene_dict[gtf_dict['gene_name']] += [i for i in range(int(fields[3])-1, int(fields[4]))]
                        
        gene_lengths = []
        for gene in gene_dict:
            if sum_type == 'gene_all':
                gene_lengths.append('{0},{1}\n'.format(gene, len(gene_dict[gene])))
            elif sum_type == 'gene_unique':
                gene_lengths.append('{0},{1}\n'.format(gene, len(set(gene_dict[gene]))))
            else:
                raise IOError('sum_type must be gene_all, gene_unique or transcript_longest')
    else:
        gene_dict = {}
        gene_count = 0
        with open(gtf_file) as input_fh:
            for line in input_fh:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'exon':
                        gtf_dict = parse_gtf_dict(fields[8])
                        if not gtf_dict['gene_name'] in gene_dict:
                            gene_count += 1
                            if gene_count % 1000 == 0:
                                print(' - {0} genes processed'.format(gene_count))
                            gene_dict[gtf_dict['gene_name']] = {}
                        if not gtf_dict['transcript_id'] in gene_dict[gtf_dict['gene_name']]:
                            gene_dict[gtf_dict['gene_name']][gtf_dict['transcript_id']] = 0
                        gene_dict[gtf_dict['gene_name']][gtf_dict['transcript_id']] += int(fields[4]) - int(fields[3]) + 1
                        
        gene_lengths = []
        for gene in gene_dict:
            transcript_lengths = gene_dict[gene].values()
            gene_lengths.append('{0},{1}\n'.format(gene, max(transcript_lengths)))

    with open(outfile, 'w') as output_fh:
        for gene in gene_lengths:
            output_fh.write(gene)


def get_transcripts(ens_gene_id, output_dir, gtf_file, fasta_str, 
    protein=False):
    exons = {}
    if protein == True:
        feature_type = 'CDS'
    else:
        feature_type = 'exon'
    with open(gtf_file) as input_fh:
        for line in input_fh:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if fields[2] == feature_type:
                    gtf_dict = parse_gtf_dict(fields[8])
                    if gtf_dict['gene_id'] == ens_gene_id:
                        if not gtf_dict['transcript_id'] in exons:
                            exons[gtf_dict['transcript_id']] = []
                        exons[gtf_dict['transcript_id']].append({
                            'chrom': fields[0],
                            'start': int(fields[3]),
                            'end': int(fields[4]),
                            'strand': fields[6],
                        })
    for transcript in exons:
        if set([i['strand'] for i in exons[transcript]]) == set('+'):
            strand = '+'
            sorted_exons = sorted(exons[transcript], key=lambda i:i['start'])
        elif set([i['strand'] for i in exons[transcript]]) == set('-'):
            strand = '-'
            sorted_exons = sorted(exons[transcript], key=lambda i:i['start'], reverse=True)
        sequence = ''
        with gzip.open(fasta_str.format(exons[transcript][0]['chrom']), 'rt') as input_fh:
            for record in SeqIO.parse(input_fh, 'fasta'):
                chrom_seq = str(record.seq)
        for exon in sorted_exons:
            if strand == '+':
                sequence += chrom_seq[exon['start']-1:exon['end']]
            elif strand == '-':
                sequence += rc(chrom_seq[exon['start']-1:exon['end']])
        if protein == True:
            sequence = translate(sequence)
        with open(os.path.join(output_dir, transcript + '.fa'), 'w') as output_fh:
            output_fh.write('>{0}\n{1}\n'.format(transcript, sequence))
        print('fasta file written: {0}/{1}.fa'.format(output_dir, transcript))


if __name__ == '__main__':
    if len(sys.argv) == 3:
        # gtf wget from ensembl (wget ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz)
        #gtf_to_gene_length(gtf_file='Mus_musculus.GRCm38.95.gtf', 
        #                   outfile='mm10_gene_lengths.txt')
        gtf_to_gene_length(gtf_file=sys.argv[1], 
                           outfile=sys.argv[2])
    elif len(sys.argv) == 4:
        gtf_to_gene_length(gtf_file=sys.argv[1], 
                           outfile=sys.argv[2],
                           include_all=sys.argv[3])
    else:
        pass
