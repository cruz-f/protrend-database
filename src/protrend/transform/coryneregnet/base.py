from typing import Dict

import pandas as pd

from protrend.io import read_from_stack, read_csv
from protrend.transform import Transformer, Connector
from protrend.utils import SetList


class CoryneRegNetTransformer(Transformer):
    default_source = 'coryneregnet'
    default_version = '0.0.0'
    default_regulation_stack = {'bsub': 'bsub_regulation.csv',
                                'cglu': 'cglu_regulation.csv',
                                'ecol': 'ecol_regulation.csv',
                                'mtub': 'mtub_regulation.csv'}

    default_regulation_columns = SetList(['TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                                          'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                                          'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source'])

    taxa_to_organism_code = {'bsub': '224308',
                             'cglu': '196627',
                             'ecol': '511145',
                             'mtub': '83332'}

    def _build_regulations(self, regulation_stack: Dict[str, str] = None) -> pd.DataFrame:

        if not regulation_stack:
            regulation_stack = self.default_regulation_stack

        dfs = []
        for file in regulation_stack:
            df = read_from_stack(stack=regulation_stack,
                                 file=file,
                                 default_columns=self.default_regulation_columns,
                                 reader=read_csv)
            df['taxonomy'] = self.taxa_to_organism_code[file]
            dfs.append(df)

        return pd.concat(dfs, axis=0)


class CoryneRegNetConnector(Connector):
    default_source = 'coryneregnet'
    default_version = '0.0.0'
