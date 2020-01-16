import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict
import time

from seqtools import alignment
from seqtools.annotation import get_sequence_from_genbank, get_sequence_from_uniprot
from seqtools.general import load_csv
from seqtools import retrieve_entrez_uniprot as reu


def _check_dir(dir_path):
    if not(os.path.exists(dir_path)):
        os.mkdir(dir_path)
    return(dir_path)


def protein(accessions_csv, out_dir, expt_name):
    aln_dir = _check_dir(os.path.join(out_dir, 'alignments'))
    ref_dir = _check_dir(os.path.join(out_dir, 'ref_sequences'))
    accessions = load_csv(accessions_csv)[1]
    sequences = []
    for gene in accessions:
        gene_name = gene['prefix'] + gene['name']
        if gene['uniprot_protein']:
            if os.path.exists(gene['uniprot_protein']):
                if gene['uniprot_protein'].endswith('.fa'):
                    fasta = SeqIO.read(gene['uniprot_protein'], 'fasta')
                    sequences.append([gene_name, str(fasta.seq)])
            else:
                transcript = reu.get_uniprot(
                    accession=gene['uniprot_protein'],
                )
                out_file = os.path.join(
                    ref_dir, '{0}_{1}.prot'.format(gene_name, gene['uniprot_protein'])
                )
                with open(out_file, 'w') as output_fh:
                    output_fh.write(transcript)
                sequences.append([gene_name, get_sequence_from_uniprot(out_file)])
                print(' - Uniprot file written: {0}'.format(out_file))
        elif gene['ncbi_protein']:
            if os.path.exists(gene['ncbi_protein']):
                if gene['ncbi_protein'].endswith('.fa'):
                    fasta = SeqIO.read(gene['ncbi_protein'], 'fasta')
                    sequences.append([gene_name, str(fasta.seq)])
            else:
                out_file = os.path.join(
                    ref_dir, '{0}_{1}.gp'.format(gene_name, gene['ncbi_protein'])
                )
                if not os.path.exists(out_file):
                    transcript = reu.get_data_entrez(
                        accession=gene['ncbi_protein'],
                        database='protein', 
                        output_format='gp'
                    )
                    time.sleep(0.2)
                    with open(out_file, 'w') as output_fh:
                        output_fh.write(transcript)
                sequences.append([gene_name, get_sequence_from_genbank(out_file)])
                print(' - GenPept file written: {0}'.format(out_file))
    if len(sequences) > 2:
        _ = alignment.run_clustalw2(
            input_seqs=sequences,
            output_file=os.path.join(aln_dir, '{0}_protein.aln'.format(expt_name))
        )
    elif len(sequences) == 2:
        _ = alignment.run_needle(
            input_seqs=sequences,
            output_file=os.path.join(aln_dir, '{0}_protein.needle'.format(expt_name))
        )


def mrna(accessions_csv, out_dir, expt_name, gene_id_name, 
         cds=False, min_contig=10,
         colour='#99ff33', conserved_all=False):
    """
    Run mRNA Orthology Pipeline:
     - ClustalO multiple alignment of all sequences with ncbi_transcript accession
     - pairwise alignment of all sequences with first sequence
    Inputs:

    Outputs:

    """
    aln_dir = _check_dir(os.path.join(out_dir, 'alignments'))
    ref_dir = _check_dir(os.path.join(out_dir, 'ref_sequences'))
    accessions = load_csv(accessions_csv)[1]
    sequences = []
    for gene in accessions:
        gene_name = gene['prefix'] + gene['name']
        if gene['ncbi_transcript']:
            if os.path.exists(gene['ncbi_transcript']):
                if gene['ncbi_transcript'].endswith('.fa'):
                    fasta = SeqIO.read(gene['ncbi_transcript'], 'fasta')
                    sequences.append([gene_name, str(fasta.seq)])
            else:
                out_file = os.path.join(
                    ref_dir, '{0}_{1}.gb'.format(gene_name, gene['ncbi_transcript'])
                )
                if not os.path.exists(out_file):
                    transcript = reu.get_data_entrez(
                                accession=gene['ncbi_transcript'],
                                database='nucleotide', 
                                output_format='gb'
                    )
                    time.sleep(0.2)
                    with open(out_file, 'w') as output_fh:
                        output_fh.write(transcript)
                sequences.append([gene_name, get_sequence_from_genbank(out_file, cds)])
                print(' - GenBank file written: {0}'.format(out_file))
    if len(sequences) > 2:
        multi_aln = alignment.run_clustalw2(
            input_seqs=sequences,
            output_file=os.path.join(aln_dir, '{0}_mRNA.aln'.format(expt_name))
        )
    contig_files = []
    for i in range(1, len(sequences)):
        aln = alignment.run_needle(
            input_seqs=[sequences[0], sequences[i]], 
            output_file=os.path.join(
                aln_dir, '{0}_{1}_mRNA.needle'.format(sequences[0][0], sequences[i][0])
            )
        )
        contig_file = os.path.join(
            aln_dir, '{0}_{1}_mRNA_contigs.csv'.format(sequences[0][0], sequences[i][0])
        )
        alignment.get_contigs_from_clustal(
            alignment=aln,
            csv_out=contig_file,
            base_sequence=sequences[0][0],
        )
        contig_files.append(contig_file)
    gene = accessions[0]
    gene_name = gene['prefix'] + gene['name']
    gb_parsed = SeqIO.read(os.path.join(ref_dir, '{0}_{1}.gb'.format(
        gene_name, gene['ncbi_transcript'],
    )), 'gb')
    if not cds:
        for contig_file in contig_files:
            headers, homo_contig = load_csv(contig_file)
            homo_name = os.path.split(contig_file)[1].split('_')[1]
            for region in homo_contig:
                if int(region['length']) > int(min_contig):
                    hc_feat = SeqFeature(
                        FeatureLocation(
                            int(region['start_wrt_{0}'.format(gene_name)])-1, 
                            int(region['end_wrt_{0}'.format(gene_name)]),
                        ), type="misc_feature",
                        qualifiers=OrderedDict([
                            ('gene', [gene_id_name]),
                            ('label', ['{0} Homology'.format(homo_name)]),
                            ('ApEinfo_revcolor', [colour]),
                            ('ApEinfo_fwdcolor', [colour]),
                        ])
                    )
                    gb_parsed.features.append(hc_feat)
            output_gb = os.path.join(out_dir, '{0}_orthology.gb'.format(expt_name))
            with open(output_gb, 'w') as output_fh:
                SeqIO.write(gb_parsed, output_fh, 'gb')


def pipeline(accessions_csv, out_dir, expt_name, gene_id_name,
             cds=False, min_contig=10, colour='#99ff33', conserved_all=False):
    protein(accessions_csv, out_dir, expt_name)
    mrna(accessions_csv, out_dir, expt_name, gene_id_name, 
         cds, min_contig, colour, conserved_all)


if __name__ == '__main__':
    pass
