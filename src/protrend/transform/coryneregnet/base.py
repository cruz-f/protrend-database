from abc import abstractmethod

import pandas as pd

from protrend.io import read_csv
from protrend.transform import MultiStackTransformer, Connector
from protrend.utils import SetList, MultiStack


class CoryneRegNetTransformer(MultiStackTransformer, source='coryneregnet', version='0.0.0', register=False):
    _taxa = ['224308',
             '196627',
             '511145',
             '83332']
    _source = ['coryneregnet'] * 4
    _net_reader = [read_csv] * 4

    default_transform_stack = {
        'network': MultiStack(
            stack=['bsub_regulation.csv',
                   'cglu_regulation.csv',
                   'ecol_regulation.csv',
                   'mtub_regulation.csv'],
            taxa=_taxa,
            source=_source,
            reader=_net_reader
        ),
    }

    default_network_columns = SetList(['TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                                       'PMID', 'Source'])

    @abstractmethod
    def transform(self):
        pass

    def merge_annotations(self, annotations: pd.DataFrame, coryneregnet: pd.DataFrame):
        df = pd.merge(annotations, coryneregnet, on='input_value', suffixes=('_annotation', '_coryneregnet'))

        # merge loci
        df = self.merge_loci(df=df, left_suffix='_annotation', right_suffix='_coryneregnet')

        # merge name
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_coryneregnet')
        return df


class CoryneRegNetConnector(Connector, source='coryneregnet', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
