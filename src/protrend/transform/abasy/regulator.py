import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Regulator
from protrend.transform.abasy.base import AbasyTransformer
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.processors import apply_processors, rstrip, lstrip
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
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'source_taxonomy'])

    def _transform_networks(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = networks.dropna(subset=['source', 'taxonomy'])
        networks['source_taxonomy'] = networks['source'] + networks['taxonomy']

        networks = self.drop_duplicates(df=networks, subset=['source_taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['source_taxonomy'])

        networks = apply_processors(networks,
                                    source_taxonomy=[rstrip, lstrip],
                                    source=[rstrip, lstrip])

        networks['mechanism'] = None
        return networks

    def transform(self):
        networks = self._build_networks()
        regulator = self._transform_networks(networks)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                                   'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                                   'sequence', 'strand', 'start', 'stop', 'gene_name_taxonomy')

        df = pd.merge(regulator, gene, left_on='source_taxonomy', right_on='gene_name_taxonomy')

        self._stack_transformed_nodes(df)

        return df
