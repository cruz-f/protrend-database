import os
from typing import Dict

import pandas as pd

from protrend.io import read_from_stack, read_csv
from protrend.transform import Transformer, Connector
from protrend.utils.processors import take_first, to_set_list
from protrend.utils import SetList, Settings


class CoryneRegNetTransformer(Transformer):
    default_source = 'coryneregnet'
    default_version = '0.0.0'
    default_regulation_stack = {'bsub': 'bsub_regulation.csv',
                                'cglu': 'cglu_regulation.csv',
                                'ecol': 'ecol_regulation.csv',
                                'mtub': 'mtub_regulation.csv'}

    default_operon_stack = {'bsub': 'bsub_operon.csv',
                            'cglu': 'cglu_operon.csv',
                            'ecol': 'ecol_operon.csv',
                            'mtub': 'mtub_operon.csv'}

    default_regulation_columns = SetList(['TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                                          'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                                          'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source'])

    default_operon_columns = SetList(['Operon', 'Orientation', 'Genes'])

    taxa_to_organism_code = {'bsub': '224308',
                             'cglu': '196627',
                             'ecol': '511145',
                             'mtub': '83332'}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._regulation_stack = {}
        self._operon_stack = {}

        self.load_regulation_stack()
        self.load_operon_stack()

    def load_regulation_stack(self, regulation_stack: Dict[str, str] = None):

        self._regulation_stack = {}

        if not regulation_stack:
            regulation_stack = self.default_regulation_stack

        for key, file in regulation_stack.items():

            sa_file = os.path.join(Settings.STAGING_AREA_PATH, self.source, self.version, file)
            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(sa_file):

                self._regulation_stack[key] = sa_file

            else:

                self._regulation_stack[key] = dl_file

    def load_operon_stack(self, operon_stack: Dict[str, str] = None):

        self._operon_stack = {}

        if not operon_stack:
            operon_stack = self.default_operon_stack

        for key, file in operon_stack.items():

            sa_file = os.path.join(Settings.STAGING_AREA_PATH, self.source, self.version, file)
            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(sa_file):

                self._operon_stack[key] = sa_file

            else:

                self._operon_stack[key] = dl_file

    @property
    def regulation_stack(self):
        return self._regulation_stack

    @property
    def operon_stack(self):
        return self._operon_stack

    def _build_regulations(self, regulation_stack: Dict[str, str] = None) -> pd.DataFrame:

        if not regulation_stack:
            regulation_stack = self.regulation_stack

        dfs = []
        for file in regulation_stack:
            df = read_from_stack(stack=regulation_stack,
                                 file=file,
                                 default_columns=self.default_regulation_columns,
                                 reader=read_csv)
            df['taxonomy'] = self.taxa_to_organism_code[file]
            dfs.append(df)

        return pd.concat(dfs, axis=0)

    def _read_operon(self, stack, file):
        fixed_names = ['Operon', 'Orientation']
        genes_names = [f'Genes{i}' for i in range(150)]
        names = fixed_names + genes_names

        df = read_from_stack(stack=stack,
                             file=file,
                             default_columns=self.default_operon_columns,
                             reader=read_csv,
                             names=names,
                             engine='python',
                             skiprows=1)
        df['Operon'] = df['Operon'].str.replace('>', '')

        df['Genes'] = df.iloc[:, 2:].values.tolist()

        df = df[['Operon', 'Orientation', 'Genes']]
        df = df.explode(column='Genes')
        df = df.dropna(subset=['Genes'])

        df = Transformer.group_by(df,
                                  column='Operon',
                                  aggregation={'Orientation': take_first},
                                  default=to_set_list)
        return df

    def _build_operons(self, operon_stack: Dict[str, str] = None) -> pd.DataFrame:

        if not operon_stack:
            operon_stack = self.operon_stack

        dfs = []
        for file in operon_stack:
            df = self._read_operon(operon_stack, file)
            df['taxonomy'] = self.taxa_to_organism_code[file]
            dfs.append(df)

        return pd.concat(dfs, axis=0)


class CoryneRegNetConnector(Connector):
    default_source = 'coryneregnet'
    default_version = '0.0.0'
