import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model import Regulator
from protrend.transform.abasy.base import AbasyTransformer, read_abasy_network
from protrend.transform.abasy.gene import GeneTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip
from protrend.utils import SetList


class RegulatorTransformer(AbasyTransformer,
                           source='abasy',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    default_transform_stack = {'gene': 'integrated_gene.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop', 'mechanism',
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'regulator_taxonomy'])

    def transform_networks(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = self.drop_duplicates(df=networks, subset=['source', 'taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['source', 'taxonomy'])
        networks = self.drop_empty_string(networks, 'source')

        networks = apply_processors(networks, source=[rstrip, lstrip])

        regulator_taxonomy = networks['source'] + networks['taxonomy']

        networks = networks.assign(regulator_taxonomy=regulator_taxonomy,
                                   mechanism=None)
        return networks

    def transform(self):
        networks = self.contact_stacks(stack=self.network_stack,
                                       taxa=self.taxa_to_organism_code,
                                       default_columns=self.default_network_columns,
                                       reader=read_abasy_network)
        regulator = self.transform_networks(networks)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                                   'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                                   'sequence', 'strand', 'start', 'stop', 'gene_taxonomy')

        df = pd.merge(regulator, gene, left_on='regulator_taxonomy', right_on='gene_taxonomy')

        df = df.dropna(subset=['locus_tag'])
        df = self.drop_empty_string(df, 'locus_tag')
        df = self.drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)

        self.stack_transformed_nodes(df)

        return df
