import os
import pytest
import shutil
from tempfile import mkdtemp

import seqtools.general
import seqtools.orthology
import seqtools.retrieve_entrez_uniprot
import seqtools.get_gene_info
import seqtools.alignment
import seqtools.annotation


def test_general(accessions_csv):
    tmp_dir = mkdtemp()
    headers, lines = seqtools.general.load_csv(
        input_csv=accessions_csv,
        sep=',',
        headers=None,
    )
    assert len(lines) == 3
    assert len(headers) == 6
    seqtools.general.write_csv(
        filename=os.path.join(tmp_dir, 'written.csv'),
        dict_list=lines,
    )
    assert os.path.exists(os.path.join(tmp_dir, 'written.csv'))
    dna1 = 'AACGAGATCCTCCATAAATAG'
    pep1 = 'NEILHK*'
    assert seqtools.general.translate(dna1) == pep1
    assert seqtools.general.transcribe(dna1) == 'AACGAGAUCCUCCAUAAAUAG'
    assert seqtools.general.rc(dna1) == 'CTATTTATGGAGGATCTCGTT'
    shutil.rmtree(tmp_dir)


def test_retrieve_entrez_uniprot(genbank_file, uniprot_file, genpept_file):
    exons = seqtools.retrieve_entrez_uniprot.get_exons_from_gb(
        gb_file=genbank_file,
    )
    assert len(exons) == 14
    exons = seqtools.retrieve_entrez_uniprot.get_exons_from_gb(
        gb_file=genbank_file,
        cds=False,
    )
    assert len(exons) == 16
    domains = seqtools.retrieve_entrez_uniprot.get_domain_details_from_gb(
        gb_file=genbank_file, 
        cds=True
    )
    assert len(domains) == 21
    features = seqtools.retrieve_entrez_uniprot.get_domain_details_from_protein(
        prot_file=uniprot_file, 
        file_type='Uniprot',
    )
    assert len(features) == 6
    features = seqtools.retrieve_entrez_uniprot.get_domain_details_from_protein(
        prot_file=genpept_file, 
        file_type='GenPept',
    )
    assert len(features) == 27
    # need to figure out how to do mock requests...


def test_get_gene_info(gtf_file, fasta_str):
    tmp_dir = mkdtemp()
    seqtools.get_gene_info.gtf_to_gene_length(
        gtf_file=gtf_file, 
        outfile=os.path.join(tmp_dir, 'gene_lengths1.csv'), 
        sum_type='transcript_longest',
    )
    _, gene_lengths1 = seqtools.general.load_csv(os.path.join(tmp_dir, 'gene_lengths1.csv'), headers=['gene', 'length'])
    assert len(gene_lengths1) == 60
    seqtools.get_gene_info.gtf_to_gene_length(
        gtf_file=gtf_file, 
        outfile=os.path.join(tmp_dir, 'gene_lengths2.csv'), 
        sum_type='gene_all',
    )
    _, gene_lengths2 = seqtools.general.load_csv(os.path.join(tmp_dir, 'gene_lengths2.csv'), headers=['gene', 'length'])
    assert len(gene_lengths2) == 60
    seqtools.get_gene_info.gtf_to_gene_length(
        gtf_file=gtf_file, 
        outfile=os.path.join(tmp_dir, 'gene_lengths3.csv'), 
        sum_type='gene_unique',
    )
    _, gene_lengths3 = seqtools.general.load_csv(os.path.join(tmp_dir, 'gene_lengths3.csv'), headers=['gene', 'length'])
    assert len(gene_lengths3) == 60
    for i in range(len(gene_lengths1)):
        assert int(gene_lengths1[i]['length']) <= int(gene_lengths3[i]['length']) <= int(gene_lengths2[i]['length'])
    nfiles1 = len(os.listdir(tmp_dir))
    seqtools.get_gene_info.get_transcripts(
        ens_gene_id='ENSG00000186092',
        output_dir=tmp_dir,
        gtf_file=gtf_file,
        fasta_str=fasta_str,
    )
    nfiles2 = len(os.listdir(tmp_dir))
    assert nfiles2 - nfiles1 == 2
    seqtools.get_gene_info.get_transcripts(
        ens_gene_id='ENSG00000186092',
        output_dir=tmp_dir,
        protein=True,
        gtf_file=gtf_file,
        fasta_str=fasta_str,
    )
    nfiles3 = len(os.listdir(tmp_dir))
    assert nfiles3 - nfiles2 == 0
    shutil.rmtree(tmp_dir)


def test_alignment(amino_acids, nucleotides):
    tmp_dir = mkdtemp()
    alignment = seqtools.alignment.run_clustalw2(
        input_seqs=amino_acids, 
        tmp_dir=tmp_dir,
    )
    assert len(alignment) == 4
    assert len(alignment[0].seq) == 81
    alignment = seqtools.alignment.run_needle(
        input_seqs=nucleotides, 
        tmp_dir=tmp_dir,
    )
    assert len(alignment) == 2
    assert len(alignment[0].seq) == 83
    seqtools.alignment.get_contigs_from_clustal(
        alignment=alignment, 
        csv_out=os.path.join(tmp_dir, 'contigs.csv'),
    )
    _, contigs = seqtools.general.load_csv(os.path.join(tmp_dir, 'contigs.csv'))
    assert len(contigs) == 9
    shutil.rmtree(tmp_dir)


def test_annotation(genbank_file, uniprot_file, sequences_file):
    tmp_dir = mkdtemp()
    seq1 = seqtools.annotation.get_sequence_from_genbank(genbank_file)
    assert len(seq1) == 3650
    seq2 = seqtools.annotation.get_sequence_from_genbank(genbank_file, cds=True)
    assert len(seq2) == 2823
    seq3 = seqtools.annotation.get_sequence_from_uniprot(uniprot_file)
    assert len(seq3) == 940
    seq4 = seqtools.annotation.clean_sequence('G;ACGmUAC*UAGC-mUAGC;UAG--CU*AGCUAmCG')
    assert seq4 == 'GACGTACTAGCTAGCTAGCTAGCTACG'
    feature = seqtools.annotation.find_sirna_and_annotate_genbank(
        ref_seq='caaggtgactgttaaatctgaaaacctcaaggttataaaggatgaagccctcagcgatggggatgacctcagggactttccaagtgacctcaagaaggcacaccatctgaagagaggggctaccatgaat',
        sequence='gactttccaagtgacctcaa',
        feature_note='siRNA1',
    )
    assert int(feature.location.start) == 73
    assert int(feature.location.end) == 93
    seqtools.annotation.annotate(
        sequence_file=sequences_file, 
        genbank_file=genbank_file, 
        output_gb=os.path.join(tmp_dir, 'annotated.gb')
    )
    domains = seqtools.retrieve_entrez_uniprot.get_domain_details_from_gb(
        gb_file=os.path.join(tmp_dir, 'annotated.gb'),
        cds=False,
    )
    assert len(domains) == 24
    shutil.rmtree(tmp_dir)
