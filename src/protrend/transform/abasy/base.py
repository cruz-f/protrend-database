import os
from typing import Dict

import pandas as pd

from protrend.io import read_json, read_from_stack, read_csv
from protrend.transform import Transformer, Connector
from protrend.utils import SetList, Settings


def read_abasy_network(file_path: str, **kwargs) -> pd.DataFrame:
    network = read_json(file_path)

    network_edges = network['elements']['edges']

    network_data = [edge['data'] for edge in network_edges]

    network_df = pd.DataFrame(network_data)
    return network_df


class AbasyTransformer(Transformer):
    default_source = 'abasy'
    default_version = '0.0.0'
    default_network_stack = {'bsub': 'bsub_network.json',
                             'cglu': 'cglu_network.json',
                             'ecol': 'ecol_network.json',
                             'mtub': 'mtub_network.json',
                             'paer': 'paer_network.json',
                             'saur': 'saur_network.json',
                             'scoe': 'scoe_network.json',
                             'spne': 'spne_network.json',
                             'spyo': 'spyo_network.json'}

    default_gene_stack = {'bsub': 'bsub_genes.tsv',
                          'cglu': 'cglu_genes.tsv',
                          'ecol': 'ecol_genes.tsv',
                          'mtub': 'mtub_genes.tsv',
                          'paer': 'paer_genes.tsv',
                          'saur': 'saur_genes.tsv',
                          'scoe': 'scoe_genes.tsv',
                          'spne': 'spne_genes.tsv',
                          'spyo': 'spyo_genes.tsv'}

    default_network_columns = SetList(['id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy'])

    default_gene_columns = SetList(['Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                                    'Product_function', 'NDA_component', 'taxonomy'])

    taxa_to_organism_code = {'bsub': '224308',
                             'cglu': '196627',
                             'ecol': '511145',
                             'mtub': '83332',
                             'paer': '208964',
                             'saur': '367830',
                             'scoe': '100226',
                             'spne': '406558',
                             'spyo': '186103'}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._network_stack = {}
        self._gene_stack = {}

        self.load_network_stack()
        self.load_gene_stack()

    def load_network_stack(self, network_stack: Dict[str, str] = None):

        self._network_stack = {}

        if not network_stack:
            network_stack = self.default_network_stack

        for key, file in network_stack.items():

            sa_file = os.path.join(Settings.STAGING_AREA_PATH, self.source, self.version, file)
            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(sa_file):

                self._network_stack[key] = sa_file

            else:

                self._network_stack[key] = dl_file

    def load_gene_stack(self, gene_stack: Dict[str, str] = None):

        self._gene_stack = {}

        if not gene_stack:
            gene_stack = self.default_gene_stack

        for key, file in gene_stack.items():

            sa_file = os.path.join(Settings.STAGING_AREA_PATH, self.source, self.version, file)
            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(sa_file):

                self._gene_stack[key] = sa_file

            else:

                self._gene_stack[key] = dl_file

    @property
    def network_stack(self):
        return self._network_stack

    @property
    def gene_stack(self):
        return self._gene_stack

    def _build_networks(self, network_stack: Dict[str, str] = None) -> pd.DataFrame:

        if not network_stack:
            network_stack = self.network_stack

        dfs = []
        for file in network_stack:
            df = read_from_stack(stack=network_stack,
                                 file=file,
                                 default_columns=self.default_network_columns,
                                 reader=read_abasy_network)
            df['taxonomy'] = self.taxa_to_organism_code[file]
            dfs.append(df)

        return pd.concat(dfs, axis=0)

    def _build_genes(self, gene_stack: Dict[str, str] = None) -> pd.DataFrame:

        if not gene_stack:
            gene_stack = self.gene_stack

        dfs = []
        for file in gene_stack:
            df = read_from_stack(stack=gene_stack,
                                 file=file,
                                 default_columns=self.default_gene_columns,
                                 reader=read_csv,
                                 sep='\t')
            df['taxonomy'] = self.taxa_to_organism_code[file]
            dfs.append(df)

        return pd.concat(dfs, axis=0)


class AbasyConnector(Connector):
    default_source = 'abasy'
    default_version = '0.0.0'
