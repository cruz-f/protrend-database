from abc import abstractmethod

import pandas as pd

from protrend.io import read_csv, read
from protrend.transform import Transformer, Connector
from protrend.transform.transformations import merge_columns, merge_loci, locus_tag_curation
from protrend.utils import apply_processors
from protrend.utils.processors import to_int_str

CORYNEREGNET_NETWORK = ['bsub_regulation.csv',
                        'cglu_regulation.csv',
                        'ecol_regulation.csv',
                        'mtub_regulation.csv']

CORYNEREGNET_TAXA = ['224308',
                     '196627',
                     '511145',
                     '83332']


def read_coryneregnet_networks(source: str, version: str) -> pd.DataFrame:
    default = pd.DataFrame(columns=['TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                                    'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                                    'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                                    'PMID', 'Source'])
    dfs = []
    for file, taxon in zip(CORYNEREGNET_NETWORK, CORYNEREGNET_TAXA):
        df = read(source=source, version=version, file=file, reader=read_csv, default=default.copy())
        df = df.assign(taxonomy=taxon, source='coryneregnet')
        dfs.append(df)

    final_df = pd.concat(dfs)
    final_df = final_df.reset_index(drop=True)
    return final_df


class CoryneRegNetTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass

    @staticmethod
    def merge_annotations(annotations: pd.DataFrame, coryneregnet: pd.DataFrame):
        df = pd.merge(annotations, coryneregnet, on='input_value', suffixes=('_annotation', '_coryneregnet'))

        # merge loci
        df = merge_loci(df=df, left_suffix='_annotation', right_suffix='_coryneregnet')

        # merge name
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_coryneregnet')

        df = apply_processors(df, taxonomy=to_int_str)
        df = locus_tag_curation(df)
        return df


class CoryneRegNetConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
