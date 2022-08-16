from abc import abstractmethod
from functools import partial

import pandas as pd

from protrend.io import read_json, read_csv, read
from protrend.transform import Transformer, Connector


ABASY_GENES = ['bsub_genes.tsv',
               'cglu_genes.tsv',
               'ecol_genes.tsv',
               'mtub_genes.tsv',
               'paer_genes.tsv',
               'saur_genes.tsv',
               'scoe_genes.tsv',
               'spne_genes.tsv',
               'spyo_genes.tsv']

ABASY_NETWORK = ['bsub_network.json',
                 'cglu_network.json',
                 'ecol_network.json',
                 'mtub_network.json',
                 'paer_network.json',
                 'saur_network.json',
                 'scoe_network.json',
                 'spne_network.json',
                 'spyo_network.json']

ABASY_TAXA = ['224308',
              '196627',
              '511145',
              '83332',
              '208964',
              '367830',
              '100226',
              '406558',
              '186103']


def _read_abasy(files, source, version, reader, default):
    dfs = []
    for file, taxon in zip(files, ABASY_TAXA):
        df = read(source=source, version=version, file=file, reader=reader, default=default.copy())
        df = df.assign(taxonomy=taxon, source='abasy')
        dfs.append(df)

    final_df = pd.concat(dfs)
    final_df = final_df.reset_index(drop=True)
    return final_df


def _read_abasy_network(file_path: str) -> pd.DataFrame:
    network = read_json(file_path)

    network_edges = network['elements']['edges']

    network_data = [edge['data'] for edge in network_edges]

    network_df = pd.DataFrame(network_data)
    network_df = network_df.rename(columns={'source': 'regulator'})

    if 'Evidence' in network_df.columns:
        network_df = network_df[network_df['Evidence'].str.lower().str.strip() == 'strong']

    return network_df


def read_abasy_networks(source: str, version: str) -> pd.DataFrame:
    default = pd.DataFrame(columns=['id', 'regulator', 'target', 'Effect', 'Evidence'])
    return _read_abasy(ABASY_NETWORK, source, version, _read_abasy_network, default)


def read_abasy_genes(source: str, version: str) -> pd.DataFrame:
    default = pd.DataFrame(columns=['Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                                    'Product_function', 'NDA_component'])
    reader = partial(read_csv, sep='\t')
    return _read_abasy(ABASY_GENES, source, version, reader, default)


class AbasyTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass


class AbasyConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
