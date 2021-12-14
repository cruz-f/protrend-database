import pandas as pd

from protrend.io import read_from_stack, read_json_frame, read_from_multi_stack
from protrend.model import Regulator
from protrend.transform.abasy.base import AbasyTransformer, read_abasy_network
from protrend.transform.abasy.gene import GeneTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip
from protrend.utils import SetList, build_stack


class RegulatorTransformer(AbasyTransformer,
                           source='abasy',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'regulator_taxonomy'])

    def transform_network(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = self.drop_duplicates(df=networks, subset=['source', 'taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['source', 'taxonomy'])
        networks = self.drop_empty_string(networks, 'source')

        networks = apply_processors(networks, source=[rstrip, lstrip])

        regulator_taxonomy = networks['source'] + networks['taxonomy']

        networks = networks.assign(regulator_taxonomy=regulator_taxonomy,
                                   mechanism=None)
        return networks

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)
        regulator = self.transform_network(network)

        gene_stack = build_stack(source=self.source, version=self.version,
                                 stack_to_load={'gene': 'integrated_gene.json'}, sa=False)

        gene = read_from_stack(stack=gene_stack,
                               key='gene',
                               columns=GeneTransformer.columns,
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
