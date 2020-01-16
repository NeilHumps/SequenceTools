from Bio import Entrez, SeqIO
import requests

API_KEY = '56b4a16903579d6184fb6fece5fae1bf2c09'


def get_data_entrez(accession, database, output_format,
                    email_address='humphryes_kirilov_neil@network.lilly.com'):
    """
    Retrieve GenBank File from Entrez ID 
    """
    Entrez.email = email_address
    try:
        handle = Entrez.efetch(db=database, 
                               rettype=output_format, 
                               id=accession.strip(),
                               api_key=API_KEY)
    except ValueError:
        raise ValueError('no Entrez results retrieved for {0}'.format(accession))
    results = handle.read()
    if results:
        return results
    else:
        raise ValueError('no Entrez results retrieved for {0}'.format(accession))


def get_uniprot(accession, 
                uniprot_url='https://www.uniprot.org/uniprot/{0}.txt'):
    """
    Retrieve Uniprot Annotation file from Uniprot ID e.g. Q15858
    """
    try:
        results = requests.get(
            uniprot_url.format(accession.strip()), allow_redirects=True
        )
    except ValueError:
        raise ValueError('no Uniprot results retrieved for {0}'.format(accession))
    if results:
        return results.content.decode("utf-8")
    else:
        raise ValueError('no Uniprot results retrieved for {0}'.format(accession))


def get_exons_from_gb(gb_file, cds=True):
    gb_parsed = SeqIO.read(gb_file, 'gb')
    if cds:
        cds_loc = [i for i in gb_parsed.features if i.type == 'CDS']
        assert len(cds_loc) == 1
        start = int(cds_loc[0].location.start)
        end = int(cds_loc[0].location.end)
    else:
        start = 0
        end = len(gb_parsed)
    exons = [i for i in gb_parsed.features if i.type == 'exon']
    exon_pos = []
    for num, exon in enumerate(exons):
        if all([int(exon.location.start) >= start, 
                int(exon.location.end) <= end]):
            exon_pos.append({
                'exon': 'exon{0}'.format(num),
                'start': int(exon.location.start),
                'end': int(exon.location.end) - 1
            })
    return(exon_pos)


def get_domain_details_from_gb(gb_file, cds=True):
    gb_parsed = SeqIO.read(gb_file, 'gb')
    if cds:
        cds_loc = [i for i in gb_parsed.features if i.type == 'CDS']
        assert len(cds_loc) == 1
        start = int(cds_loc[0].location.start)
        end = int(cds_loc[0].location.end)
    else:
        start = 0
        end = len(gb_parsed)
    exons = [i for i in gb_parsed.features if i.type == 'misc_feature']
    exon_pos = []
    for num, exon in enumerate(exons):
        if all([int(exon.location.start) >= start, 
                int(exon.location.end) <= end]):
            exon_pos.append({
                'feature': exon.qualifiers['note'][0],
                'start': int(exon.location.start),
                'end': int(exon.location.end) - 1
            })
    return(exon_pos)


def get_domain_details_from_protein(prot_file, file_type='uniprot', up_feature_types=['REGION']):
    if file_type.lower() == 'uniprot':
        prot_parsed = SeqIO.read(prot_file, 'swiss')
        feats = [i for i in prot_parsed.features if any([i.type == j for j in up_feature_types])]
        feat_pos = []
        for feat in feats:
            feat_pos.append({
                'feature': '{0}_{1}'.format(
                    feat.type, feat.qualifiers['description'].split(' ')[0]
                ),
                'start': int(feat.location.start),
                'end': int(feat.location.end) - 1
            })
        return(feat_pos)
    elif file_type.lower() == 'genpept':
        prot_parsed = SeqIO.read(prot_file, 'gb')
        regions = [i for i in prot_parsed.features if i.type == 'Region']
        sites = [i for i in prot_parsed.features if i.type == 'Site']
        feat_pos = []
        for region in regions:
            feat_pos.append({
                'feature': region.qualifiers['region_name'][0],
                'start': int(region.location.start),
                'end': int(region.location.end) - 1
            })
        for site in sites:
            feat_pos.append({
                'feature': site.qualifiers['site_type'][0],
                'start': int(site.location.start),
                'end': int(site.location.end) - 1
            })
        return(feat_pos)
    else:
        raise IOError('Invalid file_type, must be Uniprot or GenPept')

if __name__ == '__main__':
    pass
