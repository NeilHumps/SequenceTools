from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict
from seqtools.general import rc, load_csv
from seqtools.retrieve_entrez_uniprot import get_domain_details_from_protein


def get_sequence_from_genbank(gb_filename, cds=False):
    """
    Return whole sequence from genbank file and name
    Will only return CDS sequence if transcript file and cds==True
    """
    gb_parsed = SeqIO.read(gb_filename, 'gb')
    if cds:
        cds_loc = [i for i in gb_parsed.features if i.type == 'CDS']
        if len(cds_loc) == 1:
            start = int(cds_loc[0].location.start)
            end = int(cds_loc[0].location.end)
            return str(gb_parsed.seq)[start:end]
        else:
            raise IOError('CDS Sequence not available or incorrect file format')
    else:
        return str(gb_parsed.seq)


def get_sequence_from_uniprot(up_filename):
    """
    Return whole sequence from uniprot file and name
    """
    up_parsed = SeqIO.read(up_filename, 'swiss')
    return str(up_parsed.seq)


def clean_sequence(sequence, allowed='ACGU', dna=True, remove_tt=True, rev_comp=False):
    """
    Output is uppercase DNA sequence
     - only use allowed nucleotides (default=ACGT)
     - if dna == True, U -> T
     - if remove_tt == True, last two TT residues are clipped off
    """
    f_seq = ''
    for nt in sequence:
        if nt.upper() in allowed:
            if nt.upper() == 'U':
                if dna:
                    f_seq += 'T'
                else:
                    f_seq += nt.upper()
            else:
                f_seq += nt.upper()
    if all([remove_tt, f_seq[-2:] == 'TT']):
        f_seq2 = f_seq[:-2]
    else:
        f_seq2 = f_seq
    if rev_comp:
        return(rc(f_seq2))
    else:
        return(f_seq2)



def find_sirna_and_annotate_genbank(ref_seq, sequence, feature_note,
                                    feature_type="misc_feature", colour='#dc143c', allow_rc=False):
    """
    Default colour is crimson
    """
    location = ref_seq.find(sequence)
    if location == -1:
        if allow_rc:
            location = ref_seq.find(rc(sequence))
            if location == -1:
                return None
        else:
            return None
    return SeqFeature(
        FeatureLocation(
            location, 
            location + len(sequence),
        ), type=feature_type,
        qualifiers=OrderedDict([
            ('label', [feature_note]),
            ('ApEinfo_revcolor', [colour]),
            ('ApEinfo_fwdcolor', [colour]),
            ('note', [feature_note]),
        ])
    )

def add_genpept_to_genbank(gb_file, gp_file, output_gb, feature_gene, prot_file_type='genpept'):
    """
    Add all peptide domains from genpept file to genbank file containing a CDS
    """
    gb_parsed = SeqIO.read(gb_file, 'gb')   
    cds_loc = [i for i in gb_parsed.features if i.type == 'CDS']
    if len(cds_loc) == 1:
        start = int(cds_loc[0].location.start)
        end = int(cds_loc[0].location.end)
    else:
        raise IOError('CDS Sequence not available or incorrect file format')
    domains = get_domain_details_from_protein(gp_file, file_type=prot_file_type)
    feature_count = 0
    for domain in domains:
        gb_parsed.features.append(
            SeqFeature(
                FeatureLocation(
                    start + (domain['start'])*3, 
                    start + (domain['end']+1)*3,
                ), type="misc_feature",
                qualifiers=OrderedDict([
                    ('gene', [feature_gene]),
                    ('note', [domain['feature']])
                ])
        ))
        feature_count += 1
    print('{0} features added'.format(feature_count))
    with open(output_gb, 'w') as output_fh:
        SeqIO.write(gb_parsed, output_fh, 'gb')


def annotate(sequence_file, genbank_file, output_gb, 
             input_na='dna', colour='#ff8c00', allow_rc=False):
    """
    find and annotate sequences on genbank file
    sequence_file = tab-separated with 'id' and 'sequence' columns
    """
    headers, sequences = load_csv(sequence_file, sep='\t')
    gb_parsed = SeqIO.read(genbank_file, 'gb')
    ref_seq = get_sequence_from_genbank(genbank_file)
    for line in sequences:
        if input_na.upper() == 'RNA':
            dna_seq = clean_sequence(line['sequence'])
        else:
            dna_seq = line['sequence'].strip()
        seq_name = '{0}_{1}'.format(
            line['id'],
            line['sequence'],
        )
        new_feat = find_sirna_and_annotate_genbank(
            ref_seq=ref_seq, 
            sequence=dna_seq, 
            feature_note=seq_name, 
            colour=colour,
            allow_rc=allow_rc,
        )
        if new_feat:
            gb_parsed.features.append(new_feat)
        else:
            print('{0} ({1}) not found'.format(seq_name, dna_seq))
    with open(output_gb, 'w') as output_fh:
        SeqIO.write(gb_parsed, output_fh, 'gb')


if __name__ == '__main__':
    pass
