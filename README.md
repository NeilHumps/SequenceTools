# Sequence Tools Python Package

`seqtools` is a suite of python modules for sequence manipulation, including:
 - general read/write csv and biological processes such as translation and reverse complement
 - retrieval of sequence files
 - extraction of sequences and metadata
 - annotation of sequences in genbank format

## Dependencies

`seqtools` requires:

 - biopython:

`pip install biopython`

 - [CLUSTALW2](http://www.clustal.org/clustal2/)
 - [Needle](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html)

(available with `module load bioinfo`)

## `general.py`

General biology tools and file manipulation.

 - `load_csv`: Lightweight delimited file loader. If headers not provided, first line will be used as header row
 - `write_csv`: Write a list of dicts to a delimited file. If headers are not provided it will use the keys of the first dict. If the keys are not consistent throughout the list it will error.
 - `translate`: translate DNA sequence into 1 or 3 letter amino acid sequence
 - `transcribe`: transcribe DNA sequence into RNA sequence
 - `rc`: return reverse complement of DNA sequence. Ambiguous characters are tolerated.

## `retrieve_entrez_uniprot.py`

Get sequence files using NCBI Entrez API for GenBank/GenPept files or retrive Uniprot flat files.

 - `get_data_entrez`: return GenBank or GenPept files using RefSeq Accession
 - `get_uniprot`: return Uniprot flat file using accession
 - `get_exons_from_gb`: return exon positions from genbank file, including just coding exon regions if cds=True
 - `get_domain_details_from_gb`: extract feature locations from GenBank file. This returns any features labelled as 'misc_feature'
 - `get_domain_details_from_protein`: return feature positions from Uniprot or GenPept file. Uniprot features defined using up_feature_types. GenPept returns all labelled with 'Region' or 'Site'


## `get_gene_info.py`

Get transcript/gene information from GTF file

 - `gtf_to_gene_length`: Get gene lengths from a GTF file, using either the longest transcript, combined length of all unique sequence or total length of all transcripts.
 - `get_transcripts`: Get fasta files for all transcript or protein sequences for an ENS gene ID

## `alignment.py`

Run global alignment algorithms

 - `run_clustalw2`: Run CLUSTALW2 multiple sequernce alignment. Return alignment object or specify output file
 - `run_needle`: Run Needle pairwise sequernce alignment. Return alignment object or specify output file
 - `get_contigs_from_clustal`: Return regions of homology using coordinates of top sequence or base_sequence if specified

## `annotation.py`

Find and add annotations to GenBank files

 - `get_sequence_from_genbank`: Return sequence from GenBank file. (whole sequence or just CDS)
 - `get_sequence_from_uniprot`: Return sequence from GenBank file. (whole sequence or just CDS)
 - `clean_sequence`: Remove modification characters - return DNA sequence
 - `find_sirna_and_annotate_genbank`: find sequence in ref_seq and return cooordinates
 - `annotate`: Takes a list of sequences, searches for them in the GenBank file and adds annotations

## `orthology.py`

Find regions of homology between orthologues, perform alignments and annotate homology regions. Requires a manually completed accessions.csv file to specify which transcripts/proteins to use

 - `protein`: retrieves GenPept files from accessions.csv file and performs multiple sequence alignment
 - `mrna`: retrieves GenBank files from accessions.csv, performs multiple sequence alignment and pairwise alignment, then annotates regions of homology on genbank file.
