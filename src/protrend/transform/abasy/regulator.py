import pandas as pd

from protrend.io.utils import read_gene
from protrend.model import Regulator
from protrend.report import ProtrendReporter
from protrend.transform.abasy.base import AbasyTransformer, read_abasy_networks
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList
from protrend.utils.constants import UNKNOWN
from protrend.utils.processors import apply_processors, rstrip, lstrip


class RegulatorTransformer(GeneMixIn, AbasyTransformer,
                           source='abasy',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop', 'mechanism',
                       'id', 'regulator', 'target', 'Effect', 'Evidence', 'source', 'taxonomy', 'regulator_taxonomy'])

    @staticmethod
    def transform_network(networks: pd.DataFrame) -> pd.DataFrame:
        networks = drop_duplicates(df=networks, subset=['regulator', 'taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['regulator', 'taxonomy'])
        networks = drop_empty_string(networks, 'regulator')

        networks = apply_processors(networks, regulator=[rstrip, lstrip])

        regulator_taxonomy = networks['regulator'] + networks['taxonomy']

        networks = networks.assign(regulator_taxonomy=regulator_taxonomy,
                                   mechanism=UNKNOWN)
        return networks

    @staticmethod
    def transform_gene(gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                              'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                              'protein_sequence', 'strand', 'start', 'stop', 'gene_taxonomy')
        return gene

    def transform(self):
        network = read_abasy_networks(self.source, self.version)
        regulator = self.transform_network(network)

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        gene = self.transform_gene(gene)

        df = pd.merge(regulator, gene, left_on='regulator_taxonomy', right_on='gene_taxonomy')

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=df.shape[0], properties=df.shape[1])

        self.stack_transformed_nodes(df)

        return df
