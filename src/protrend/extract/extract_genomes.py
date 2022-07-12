import gzip
import os
import shutil
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from protrend.bioapis import fetch_sequences
from protrend.io import write_json_frame
from protrend.log import ProtrendLogger
from protrend.utils import Settings


def _get_uniprot_accession(db_xref: List[str]) -> Optional[str]:
    """
    Retrieves the uniprot accession for the given db xref.
    """
    for xref in db_xref:
        if not xref:
            continue

        if xref.startswith('UniProtKB:'):
            return xref.split(':')[1]
    return


def read_gb_file(gb_file_path: Path):
    """
    Reads the given genbank file.
    """
    with gzip.open(gb_file_path, "rt") as handle:
        try:
            return [record for record in SeqIO.parse(handle, "genbank")]
        except Exception as e:
            ProtrendLogger.log.info(f'Error reading {gb_file_path}: {e}')
            return []


def read_genomes_ftps(genomes_file: Path = None) -> pd.DataFrame:
    """
    Reads the genomes file and returns a pandas dataframe.
    """
    if not genomes_file:
        genomes_file = Settings.genomes_ftps

    return pd.read_csv(genomes_file)


def fetch_genomes(genomes_ftps: pd.DataFrame):
    """
    Fetches the genomes from the given genomes ftps.
    """
    genomes_ftps = genomes_ftps.dropna(subset=['ncbi_taxonomy', 'genbank_ftp'])
    ncbi_taxa = genomes_ftps['ncbi_taxonomy'].to_list()
    ftp_paths = genomes_ftps['genbank_ftp'].to_list()

    output_sequences = []
    for ftp_path, ncbi_taxonomy in tqdm(zip(ftp_paths, ncbi_taxa), desc='fetch genomes from NCBI FTP',
                                        total=len(ftp_paths)):
        file_path = Settings.genomes_database.joinpath(f'{ncbi_taxonomy}.gbff.gz')
        if not file_path.exists():
            fetch_sequences(ftp_path, file_path)
        output_sequences.append((file_path, ncbi_taxonomy))

    return output_sequences


def _parse_record(record: SeqIO.SeqRecord, promoter_region_length: int) -> dict:
    """
    Parses the given record.
    """
    genome = {'locus_tag': [],
              'name': [],
              'synonyms': [],
              'uniprot_accession': [],
              'genbank_accession': [],
              'gene_sequence': [],
              'gene_start': [],
              'gene_end': [],
              'gene_strand': [],
              'protein_sequence': [],
              'promoter_sequence': [],
              'promoter_start': [],
              'promoter_end': [],
              'promoter_strand': []}

    for feature in record.features:

        if feature.type == 'CDS':
            # Identifies the strand of the gene
            strand = feature.strand

            # Identifies the start position of the gene on the sense strand (5'-3')
            gene_start = feature.location.start.position

            # Identifies the end position of the gene on the anti-sense strand (3'-5')
            gene_end = feature.location.end.position

            if strand == 1:
                promoter_region_start = int(gene_start - promoter_region_length)
                promoter_region_end = int(gene_start)

                if promoter_region_start < 0:
                    promoter_region_start = 0

                promoter_region_sequence = record.seq[promoter_region_start:promoter_region_end]

            elif strand == -1:
                promoter_region_start = int(gene_end)
                promoter_region_end = int(gene_end + promoter_region_length)

                if promoter_region_end > len(record.seq):
                    promoter_region_end = len(record.seq)

                promoter_region_sequence = record.seq[promoter_region_start:promoter_region_end].reverse_complement()

            else:
                continue

            locus_tag = feature.qualifiers.get('locus_tag', [None])[0]
            name = feature.qualifiers.get('gene', [None])[0]
            synonyms = feature.qualifiers.get('gene_synonym', [None])[0]
            uniprot_accession = _get_uniprot_accession(feature.qualifiers.get('db_xref', [None]))
            genbank_accession = feature.qualifiers.get('protein_id', [None])[0]
            gene_sequence = str(feature.extract(record.seq))
            protein_sequence = feature.qualifiers.get('translation', [None])[0]
            promoter_sequence = str(promoter_region_sequence)
            promoter_start = promoter_region_start
            promoter_end = promoter_region_end
            promoter_strand = strand

            genome['locus_tag'].append(locus_tag)
            genome['name'].append(name)
            genome['synonyms'].append(synonyms)
            genome['uniprot_accession'].append(uniprot_accession)
            genome['genbank_accession'].append(genbank_accession)
            genome['gene_sequence'].append(gene_sequence)
            genome['gene_start'].append(gene_start)
            genome['gene_end'].append(gene_end)
            genome['gene_strand'].append(strand)
            genome['protein_sequence'].append(protein_sequence)
            genome['promoter_sequence'].append(promoter_sequence)
            genome['promoter_start'].append(promoter_start)
            genome['promoter_end'].append(promoter_end)
            genome['promoter_strand'].append(promoter_strand)

    return genome


def build_genome_database(gb_file: Path,
                          promoter_region_length: int = 150) -> pd.DataFrame:
    """
    Retrieves the promoter sequences for the given organisms.
    """
    records = read_gb_file(gb_file)
    genomes = []
    for record in records:
        genome_record = _parse_record(record, promoter_region_length)
        genome_record = pd.DataFrame(genome_record)
        genomes.append(genome_record)

    if genomes:
        return pd.concat(genomes, ignore_index=True)

    return pd.DataFrame()


def build_genomes_database(genomes_file: Path = None,
                           promoter_region_length: int = 150):
    """
    Builds the genomes' database.
    """
    if not genomes_file:
        genomes_file = Settings.genomes_ftps

    genomes_ftps = read_genomes_ftps(genomes_file)

    genomes_gnbs = fetch_genomes(genomes_ftps)

    for genome_gnb, ncbi_taxonomy in tqdm(genomes_gnbs, desc='build genomes from GenBanks',
                                          total=len(genomes_gnbs)):
        genome_path = Path(Settings.genomes_database).joinpath(f'{ncbi_taxonomy}.json')
        if not genome_path.exists():
            genome = build_genome_database(genome_gnb, promoter_region_length)
            write_json_frame(genome_path, genome)


if __name__ == '__main__':
    from protrend.runners.runners import run_logger
    run_logger('extract_genomes')
    build_genomes_database()
