import os
from typing import Dict

import pandas as pd

from protrend.io import read_from_stack
from protrend.transform import Transformer, Connector
from protrend.transform.processors import apply_processors, to_set_list
from protrend.utils import SetList, Settings
from protrend.utils.miscellaneous import is_null


def read_ecol_fang_et_al_2017(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_excel(file_path, sheet_name='TU', **kwargs)


def read_mtub_turkarslan_et_al_2015(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_excel(file_path, **kwargs)


def read_paer_vasquez_et_al_2011(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_excel(file_path, skiprows=3, **kwargs)


def read_bsub_faria_et_al_2017(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_excel(file_path, sheet_name='S1 Ful Net V1', **kwargs)


def read_regulators_bsub_faria_et_al_2017(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_excel(file_path, sheet_name='S2 Regulators', skiprows=7, **kwargs)


class LiteratureTransformer(Transformer):
    default_source = 'literature'
    default_version = '0.0.0'
    default_network_stack = {'bsub': 'faria_2016.xlsx',
                             'ecol': 'fang_2017.xlsx',
                             'mtub': 'turkarslan_2015.xls',
                             'paer': 'vasquez_2011.xls'}

    taxa_to_organism_code = {'bsub': '224308',
                             'ecol': '511145',
                             'mtub': '83332',
                             'paer_PAO1': '208964',
                             'paer_PA103': '1081927',
                             'paer_PA14': '652611',
                             'paer_PAK': '1009714'}

    taxa_to_source_code = {'bsub': 'bsub_faria_et_al_2017',
                           'ecol': 'ecol_fang_et_al_2017',
                           'mtub': 'mtub_turkarslan_et_al_2015',
                           'paer': 'paer_vasquez_et_al_2011',
                           'paer_PAO1': 'paer_vasquez_et_al_2011',
                           'paer_PA103': 'paer_vasquez_et_al_2011',
                           'paer_PA14': 'paer_vasquez_et_al_2011',
                           'paer_PAK': 'paer_vasquez_et_al_2011'}

    default_network_columns = SetList(['regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
                                       'publication', 'taxonomy', 'source',
                                       'network_id'])

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._network_stack = {}

        self.load_network_stack()

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

    @property
    def network_stack(self):
        return self._network_stack

    def _build_ecol_fang_et_al_2017(self) -> pd.DataFrame:
        df = read_from_stack(stack=self.network_stack,
                             file='ecol',
                             default_columns=self.default_network_columns,
                             reader=read_ecol_fang_et_al_2017)
        df['taxonomy'] = self.taxa_to_organism_code['ecol']
        df['source'] = self.taxa_to_source_code['ecol']

        if df.empty:
            return df

        df = df.dropna(subset=['TF_id', 'gene_ids'])

        def split_regs(item):
            if is_null(item):
                return []
            return item.split(';')

        def split_genes(item):
            if is_null(item):
                return []
            return item.split(',')

        df = apply_processors(df=df, TF_id=split_regs, gene_ids=split_genes)
        df = df.explode(column='TF_id')

        df = df.drop(columns=['evidence', 'ev_level'])
        columns = {'TF': 'regulator_name',
                   'TF_id': 'regulator_locus_tag',
                   'TU_name': 'operon',
                   'genes': 'genes_name',
                   'gene_ids': 'genes_locus_tag',
                   'effect': 'regulatory_effect'}
        df = df.rename(columns=columns)

        df['mechanism'] = None
        df['evidence'] = None
        df['effector'] = None
        df['publication'] = None

        return df

    def _build_mtub_turkarslan_et_al_2015(self) -> pd.DataFrame:
        df = read_from_stack(stack=self.network_stack,
                             file='mtub',
                             default_columns=self.default_network_columns,
                             reader=read_mtub_turkarslan_et_al_2015)
        df['taxonomy'] = self.taxa_to_organism_code['mtub']
        df['source'] = self.taxa_to_source_code['mtub']

        if df.empty:
            return df

        df = self.select_columns(df, 'Transcription Factor', 'Target Gene', 'Operon', 'Differential Expression')
        df = df.dropna(subset=['Transcription Factor', 'Target Gene', 'Operon', 'Differential Expression'])

        repeated_header_mask = df['Transcription Factor'] != 'Transcription Factor'
        df = df[repeated_header_mask]

        de_mask = df['Differential Expression'] != 'no'
        df = df[de_mask]

        df_gb = Transformer.group_by(df, column='Operon', aggregation={}, default=to_set_list)
        df_gb = df_gb.drop(columns=['Transcription Factor'])

        df = df.drop(columns=['Target Gene', 'Differential Expression'])

        df = pd.merge(df, df_gb, on='Operon')

        columns = {'Transcription Factor': 'regulator_locus_tag',
                   'Operon': 'operon',
                   'Target Gene': 'genes_locus_tag',
                   'Differential Expression': 'regulatory_effect'}
        df = df.rename(columns=columns)

        df['regulator_name'] = None
        df['genes_name'] = None
        df['mechanism'] = None
        df['evidence'] = None
        df['effector'] = None
        df['publication'] = None

        return df

    def _build_paer_vasquez_et_al_2011(self) -> pd.DataFrame:
        df = read_from_stack(stack=self.network_stack,
                             file='paer',
                             default_columns=self.default_network_columns,
                             reader=read_paer_vasquez_et_al_2011)
        df['source'] = self.taxa_to_source_code['paer']

        if df.empty:
            return df

        df = df.dropna(subset=['Regulator (TF or sigma)', 'Target gene'])
        df['Operon'] = df['Operon'].fillna(df['Target gene'])

        df_gb = self.group_by(df, column='Operon', aggregation={}, default=to_set_list)
        df_gb = df_gb.drop(columns=['Regulator (TF or sigma)'])

        df = self.select_columns(df, 'Regulator (TF or sigma)', 'Operon')

        df = pd.merge(df, df_gb, on='Operon')

        df = df.explode(column='P. aeruginosa Strain')
        df = df.dropna(subset=['P. aeruginosa Strain'])

        atcc_mask = df['P. aeruginosa Strain'] != 'ATCC'
        df = df[atcc_mask]

        def split_strains(item: str) -> pd.DataFrame:
            if item == 'PA103, PA14 ans PAO1':
                return ['PA103', 'PA14', 'PAO1']

            if item == 'PAK and PAO1':
                return ['PAK', 'PAO1']

            return [item]

        df = apply_processors(df, **{'P. aeruginosa Strain': split_strains})
        df = df.explode(column='P. aeruginosa Strain')

        def set_taxonomy_by_strain(item: str) -> pd.DataFrame:
            if item == 'PAO1':
                return '208964'

            if item == 'PA103':
                return '1081927'

            if item == 'PA14':
                return '652611'

            if item == 'PAK':
                return '1009714'

        df = apply_processors(df, **{'P. aeruginosa Strain': set_taxonomy_by_strain})

        columns = {'Regulator (TF or sigma)': 'regulator_name',
                   'Operon': 'operon',
                   'Target gene': 'genes_name',
                   'mode of regulation': 'regulatory_effect',
                   'Experimental Evidence': 'evidence',
                   'PubMed Referencea': 'publication',
                   'P. aeruginosa Strain': 'taxonomy'}
        df = df.rename(columns=columns)

        df['regulator_locus_tag'] = None
        df['genes_locus_tag'] = None
        df['mechanism'] = None

        return df

    def _build_bsub_faria_et_al_2017(self) -> pd.DataFrame:
        df = read_from_stack(stack=self.network_stack,
                             file='bsub',
                             default_columns=self.default_network_columns,
                             reader=read_bsub_faria_et_al_2017)
        df['taxonomy'] = self.taxa_to_organism_code['bsub']
        df['source'] = self.taxa_to_source_code['bsub']

        if df.empty:
            return df

        df.columns = [col.rstrip().lstrip() for col in df.columns]

        df = self.select_columns(df, 'BSU Number', 'Gene_Name', 'Operon',
                                 'Sigma factor number', 'Regulator number',
                                 'Regulation sign', 'Metabolite(s) number')
        df = df.dropna(subset=['BSU Number'])

        regulator_mask = df['Sigma factor number'].notnull() | df['Regulator number'].notnull()
        df = df[regulator_mask]

        def split_regulator(item):
            if is_null(item):
                return []

            return item.split('|')

        df = apply_processors(df=df, **{'Sigma factor number': split_regulator, 'Regulator number': split_regulator})
        df = df.explode(column='Sigma factor number')
        df = df.explode(column='Regulator number')

        df['Operon'] = df['Operon'].fillna(df['Gene_Name'])

        df_gb = self.group_by(df, column='Operon', aggregation={}, default=to_set_list)
        df_gb = self.select_columns(df_gb, 'BSU Number', 'Gene_Name', 'Operon')

        df = self.select_columns(df,
                                 'Sigma factor number', 'Regulator number', 'Operon',
                                 'Regulation sign', 'Metabolite(s) number')

        df = pd.merge(df, df_gb, on='Operon')

        df_sigma = df.drop(columns=['Regulator number'])
        df_sigma = df_sigma.dropna(subset=['Sigma factor number'])
        df_sigma_mask = df_sigma['Sigma factor number'] != ''
        df_sigma = df_sigma[df_sigma_mask]
        df_sigma = df_sigma.rename(columns={'Sigma factor number': 'Regulator number'})
        df_sigma['mechanism'] = 'Sigma factor'

        df_regulator = df.drop(columns=['Sigma factor number'])
        df_regulator = df_regulator.dropna(subset=['Regulator number'])
        df_regulator_mask = df_regulator['Regulator number'] != ''
        df_regulator = df_regulator[df_regulator_mask]
        df_regulator['mechanism'] = 'Transcription factor'

        df = pd.concat([df_regulator, df_sigma], axis=0)

        df_regulators = read_from_stack(stack=self.network_stack,
                                        file='bsub',
                                        default_columns=self.default_network_columns,
                                        reader=read_regulators_bsub_faria_et_al_2017)

        df_regulators.columns = [col.rstrip().lstrip() for col in df_regulators.columns]
        df_regulators = self.select_columns(df_regulators, 'Regulator name', 'BSU', 'Number', 'Mechanism')

        # 'BSU Number', 'Gene_Name', 'Operon', 'Regulator number', 'Regulation sign', 'Involved Metabolite(s)'

        df = pd.merge(df, df_regulators, left_on='Regulator number', right_on='Number')
        df = df.drop(columns=['Regulator number', 'Number'])

        columns = {'Regulator name': 'regulator_name',
                   'BSU': 'regulator_locus_tag',
                   'Operon': 'operon',
                   'BSU Gene_Name': 'genes_name',
                   'BSU Number': 'genes_locus_tag',
                   'Regulation sign': 'regulatory_effect',
                   'Mechanism': 'mechanism',
                   'Involved Metabolite(s)': 'effector'}
        df = df.rename(columns=columns)

        df['evidence'] = None
        df['publication'] = None

        return df

    def _build_network(self) -> pd.DataFrame:
        dfs = [self._build_ecol_fang_et_al_2017(), self._build_mtub_turkarslan_et_al_2015(),
               self._build_paer_vasquez_et_al_2011(), self._build_bsub_faria_et_al_2017(), ]
        df = pd.concat(dfs, axis=0)

        df['network_id'] = df['regulator_locus_tag'] + df['regulator_name'] + df['operon'] + df['taxonomy'] + df['source']

        return df


class LiteratureConnector(Connector):
    default_source = 'literature'
    default_version = '0.0.0'
