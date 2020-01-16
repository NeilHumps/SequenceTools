import os
import pytest
from tempfile import mkdtemp


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], '../data')
WORK_DIR = mkdtemp()


@pytest.fixture(scope="module")
def accessions_csv():
    return(os.path.join(TEST_DATA_DIR, 'accessions.csv'))


@pytest.fixture(scope="module")
def genbank_file():
    return(os.path.join(TEST_DATA_DIR, 'hXPC_NM_004628.5.gb'))


@pytest.fixture(scope="module")
def genpept_file():
    return(os.path.join(TEST_DATA_DIR, 'hXPC_NP_004619.3.gp'))


@pytest.fixture(scope="module")
def uniprot_file():
    return(os.path.join(TEST_DATA_DIR, 'hXPC_Q01831.txt'))


@pytest.fixture(scope="module")
def gtf_file():
    return(os.path.join(TEST_DATA_DIR, 'human38_test.gtf'))


@pytest.fixture(scope="module")
def fasta_str():
    return(os.path.join(TEST_DATA_DIR, 'Homo_sapiens.GRCh38.dna.chromosome.{0}.fa.gz'))


@pytest.fixture(scope="module")
def amino_acids():
    return([
        ['pep1', 'lvhifllilralqlltrlvlslqpiplksatakgkkpskerltadpggssetssqvlenhtkpktskgtkqeetfa'],
        ['pep2', 'lvhifflralqlltrlvlslqpiplksatagkkpskerltadpggssetssqlenhtkktsykgtkqeetfa'],
        ['pep3', 'lvhifllralqlltrlvlslqpipalksatakgkkpskerltadpgdfhgssetssqvlenhtkpktsqeetfa'],
        ['pep4', 'lvhiflliralqlltrlvlslqpikplksatkpkerltakkdpggssetvssqvlenhtkpktskgtkqeetfa'],
    ])


@pytest.fixture(scope="module")
def nucleotides():
    return([
        ['dna1', 'TAGCTAGCTAGCTATATTAGCGCGCGTATGACTGAGAGAGATTCTCAGTAGCTAGCTCGGACGATCATGCGCAT'],
        ['dna2', 'TAGCTACAGCTAGACTGAGCGCGCGTATGACTGAGAGAGATTTCGCTACGCTCGGTGGTGACGATCCTATGCGCAT'],
    ])


@pytest.fixture(scope="module")
def sequences_file():
    return(os.path.join(TEST_DATA_DIR, 'sequences.txt'))

