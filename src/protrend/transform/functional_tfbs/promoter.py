import os
from pathlib import Path
import gzip

import pandas as pd
from Bio import SeqIO
from neo4j.exceptions import Neo4jError, DriverError

from protrend.bioapis import fetch_sequences_to_ncbi_ftp
from protrend.model import Gene, Organism, PromoterRegion
from protrend.transform.functional_tfbs.base import FunctionalTFBSTransformer


def read_gb_file(gb_file_path: Path):
    with gzip.open(gb_file_path, "rt") as handle:
        for record in SeqIO.parse(handle, 'gb'):
            return record


class PromoterRegionTransformer(FunctionalTFBSTransformer,
                                source='functional_tfbs',
                                version='0.0.0',
                                node=PromoterRegion,
                                order=100,
                                register=True):

    @staticmethod
    def fetch_genes():
        try:
            return Gene.node_to_df()
        except (Neo4jError, DriverError) as e:
            print(e)
            return pd.DataFrame()

    @staticmethod
    def fetch_organisms():
        try:
            return Organism.node_to_df()
        except (Neo4jError, DriverError) as e:
            print(e)
            return pd.DataFrame()

    def fetch_sequences(self, organisms: pd.DataFrame):
        """
        Fetches the sequences for the given organisms.
        """
        output_folder = os.path.join(self.write_path, 'sequences')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        ftp_paths = organisms['genbank_ftp'].to_list()
        protrend_ids = organisms['protrend_id'].to_list()
        return fetch_sequences_to_ncbi_ftp(protrend_ids=protrend_ids, ftp_paths=ftp_paths, output_folder=output_folder)

    @staticmethod
    def retrieve_promoter_sequences(gb_file: Path,
                                    promoter_region_length: int = 150) -> pd.DataFrame:
        """
        Retrieves the promoter sequences for the given organisms.
        :param gb_file: Path to the genbank file.
        :param promoter_region_length: Length of the promoter region.
        :return: Dataframe with the promoter sequences.
        """
        record = read_gb_file(gb_file)

        promoters = {'locus_tag': [],
                     'uniprot_accession': [],
                     'promoter_sequence': [],
                     'start': [],
                     'end': [],
                     'strand': []}
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

                    promoter_region = record.seq[promoter_region_start:promoter_region_end]

                elif strand == -1:
                    promoter_region_start = int(gene_end)
                    promoter_region_end = int(gene_end + promoter_region_length)

                    if promoter_region_end > len(record.seq):
                        promoter_region_end = len(record.seq)

                    promoter_region = record.seq[promoter_region_start:promoter_region_end].reverse_complement()

                else:
                    continue

                locus_tag = feature.qualifiers.get('locus_tag', [None])[0]
                uniprot_id = feature.qualifiers.get('protein_id', [None])[0]

                promoters['locus_tag'].append(locus_tag)
                promoters['uniprot_accession'].append(uniprot_id)
                promoters['promoter_sequence'].append(str(promoter_region))
                promoters['start'].append(promoter_region_start)
                promoters['end'].append(promoter_region_end)
                promoters['strand'].append(strand)

        return pd.DataFrame(promoters)

    def transform(self) -> pd.DataFrame:
        organisms = self.fetch_organisms()
        sequences = self.fetch_sequences(organisms)

        promoters = []
        for sequence in sequences:
            promoter = self.retrieve_promoter_sequences(sequence)
            promoters.append(promoter)

        promoters = pd.concat(promoters)
        promoters['locus_tag_merge'] = promoters['locus_tag'].str.lower()

        genes = self.fetch_genes()
        genes = genes[['locus_tag', 'uniprot_accession']]
        genes['locus_tag_merge'] = genes['locus_tag_merge'].str.lower()

        promoters = pd.merge(promoters, genes, on='locus_tag_merge', suffixes=('_promoters', '_genes'))
        promoters = promoters.drop(columns=['locus_tag_promoters', 'uniprot_accession_promoters', 'locus_tag_merge'])
        promoters = promoters.rename(columns={'locus_tag_genes': 'locus_tag',
                                              'uniprot_accession_genes': 'uniprot_accession'})

        self.stack_transformed_nodes(promoters)
        return promoters
