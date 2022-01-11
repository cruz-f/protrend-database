import pandas as pd

from protrend.io import read_from_stack, read_json_frame, read_from_multi_stack
from protrend.model import Regulator
from protrend.transform.abasy.base import AbasyTransformer
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList, build_stack
from protrend.utils.processors import apply_processors, rstrip, lstrip


class RegulatorTransformer(GeneMixIn, AbasyTransformer,
                           source='abasy',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'id', 'regulator', 'target', 'Effect', 'Evidence', 'source', 'taxonomy', 'regulator_taxonomy'])

    @staticmethod
    def transform_network(networks: pd.DataFrame) -> pd.DataFrame:
        networks = drop_duplicates(df=networks, subset=['regulator', 'taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['regulator', 'taxonomy'])
        networks = drop_empty_string(networks, 'regulator')

        networks = apply_processors(networks, regulator=[rstrip, lstrip])

        regulator_taxonomy = networks['regulator'] + networks['taxonomy']

        networks = networks.assign(regulator_taxonomy=regulator_taxonomy,
                                   mechanism=None)
        return networks

    @staticmethod
    def transform_gene(gene_stack: dict) -> pd.DataFrame:
        gene = read_from_stack(stack=gene_stack,
                               key='gene',
                               columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                              'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                              'sequence', 'strand', 'start', 'stop', 'gene_taxonomy')
        return gene

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)
        regulator = self.transform_network(network)

        gene_stack = build_stack(source=self.source, version=self.version,
                                 stack_to_load={'gene': 'integrated_gene.json'}, sa=False)
        gene = self.transform_gene(gene_stack)

        df = pd.merge(regulator, gene, left_on='regulator_taxonomy', right_on='gene_taxonomy')

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)

        self.stack_transformed_nodes(df)

        return df
