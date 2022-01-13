from abc import abstractmethod

import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.transform import MultiStackTransformer, Connector
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates, merge_columns
from protrend.utils import SetList, is_null, MultiStack
from protrend.utils.processors import apply_processors, to_int_str, rstrip, lstrip, to_set_list


def read_bsub_faria_et_al_2017(file_path: str, **kwargs) -> pd.DataFrame:
    # network
    df = pd.read_excel(file_path, sheet_name='S1 Ful Net V1', **kwargs)
    df.columns = [col.rstrip().lstrip() for col in df.columns]

    df_sigma = select_columns(df, 'BSU Number', 'Sigma factor number',
                              'Regulation sign', 'Involved Metabolite(s)')
    df_sigma = df_sigma.rename(columns={'Sigma factor number': 'Regulator number'})
    df_sigma = df_sigma.reset_index(drop=True)

    df_regulator = select_columns(df, 'BSU Number', 'Regulator number',
                                  'Regulation sign', 'Involved Metabolite(s)')
    df_regulator = df_regulator.reset_index(drop=True)

    df = pd.concat([df_regulator, df_sigma])

    def split_regulator(item):
        if is_null(item):
            return []

        return item.split('|')

    df = apply_processors(df=df, **{'Regulator number': split_regulator})
    df = df.explode(column='Regulator number')

    df = apply_processors(df, **{'Regulator number': [rstrip, lstrip, to_int_str]})

    df = df.dropna(subset=['BSU Number', 'Regulator number', 'Regulation sign'])
    df = drop_empty_string(df, 'BSU Number', 'Regulator number', 'Regulation sign')
    df = drop_duplicates(df=df,
                         subset=['BSU Number', 'Regulator number',
                                 'Regulation sign', 'Involved Metabolite(s)'],
                         perfect_match=True)

    # regulators page with BSU locus_tag
    regs = pd.read_excel(file_path, sheet_name='S2 Regulators', skiprows=7, **kwargs)
    regs.columns = [col.rstrip().lstrip() for col in regs.columns]

    regs = select_columns(regs, 'BSU', 'Number', 'Mechanism')
    regs = apply_processors(regs, **{'Number': to_int_str})

    # 'BSU Number', 'Regulator number', 'Regulation sign', 'Involved Metabolite(s)' + 'BSU', 'Number', 'Mechanism'
    df = pd.merge(df, regs, left_on='Regulator number', right_on='Number')
    df = df.drop(columns=['Regulator number', 'Number'])

    columns = {'BSU': 'regulator_locus_tag',
               'BSU Number': 'gene_locus_tag',
               'Regulation sign': 'regulatory_effect',
               'Mechanism': 'mechanism',
               'Involved Metabolite(s)': 'effector_name'}
    df = df.rename(columns=columns)
    df = df.assign(evidence=None, publication=None)
    return df


def read_ecol_fang_et_al_2017(file_path: str, **kwargs) -> pd.DataFrame:
    df = pd.read_excel(file_path, sheet_name='TU', **kwargs)

    df = apply_processors(df=df, TF_id=[lstrip, rstrip], gene_ids=[lstrip, rstrip])
    df = df.dropna(subset=['TF_id', 'gene_ids', 'effect'])
    df = drop_empty_string(df, 'TF_id', 'gene_ids', 'effect')
    df = drop_duplicates(df, subset=['TF_id', 'gene_ids', 'effect'], perfect_match=True)

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
    df = df.explode(column='gene_ids')

    df = df.dropna(subset=['TF_id', 'gene_ids', 'effect'])
    df = drop_empty_string(df, 'TF_id', 'gene_ids', 'effect')
    df = drop_duplicates(df, subset=['TF_id', 'gene_ids', 'effect'], perfect_match=True)

    df = select_columns(df, 'TF_id', 'gene_ids', 'effect')
    columns = {'TF_id': 'regulator_locus_tag',
               'gene_ids': 'gene_locus_tag',
               'effect': 'regulatory_effect'}
    df = df.rename(columns=columns)
    df = df.assign(evidence=None, effector_name=None, mechanism=None, publication=None)
    return df


def read_mtub_turkarslan_et_al_2015(file_path: str, **kwargs) -> pd.DataFrame:
    df = pd.read_excel(file_path, **kwargs)

    df = select_columns(df, 'Transcription Factor', 'Target Gene', 'Differential Expression')

    repeated_header_mask = df['Transcription Factor'] != 'Transcription Factor'
    df = df[repeated_header_mask]

    de_mask = df['Differential Expression'] != 'no'
    df = df[de_mask]

    df = df.dropna(subset=['Transcription Factor', 'Target Gene', 'Differential Expression'])
    df = drop_empty_string(df, 'Transcription Factor', 'Target Gene', 'Differential Expression')
    df = drop_duplicates(df, subset=['Transcription Factor', 'Target Gene', 'Differential Expression'],
                         perfect_match=True)

    columns = {'Transcription Factor': 'regulator_locus_tag',
               'Target Gene': 'gene_locus_tag',
               'Differential Expression': 'regulatory_effect'}
    df = df.rename(columns=columns)
    df = df.assign(evidence=None, effector_name=None, mechanism=None, publication=None)
    return df


def read_paer_vasquez_et_al_2011(file_path: str, **kwargs) -> pd.DataFrame:
    df = pd.read_excel(file_path, skiprows=3, **kwargs)

    df = select_columns(df, 'Regulator (TF or sigma)', 'Target Gene', 'mode of regulation',
                        'Experimental Evidence', 'PubMed Referencea', 'P. aeruginosa Strain')

    df = df.dropna(subset=['Regulator (TF or sigma)', 'Target Gene', 'mode of regulation', 'P. aeruginosa Strain'])
    df = drop_empty_string(df, 'Regulator (TF or sigma)', 'Target Gene', 'mode of regulation',
                           'P. aeruginosa Strain')
    df = drop_duplicates(df=df,
                         subset=['Regulator (TF or sigma)', 'Target Gene', 'mode of regulation',
                                 'P. aeruginosa Strain'],
                         perfect_match=True)

    df = df.explode(column='P. aeruginosa Strain')

    atcc_mask = df['P. aeruginosa Strain'] != 'ATCC'
    df = df[atcc_mask]

    def split_strains(item: str) -> list:
        if item == 'PA103, PA14 ans PAO1':
            return ['PA103', 'PA14', 'PAO1']

        if item == 'PAK and PAO1':
            return ['PAK', 'PAO1']

        return [item]

    df = apply_processors(df, **{'P. aeruginosa Strain': split_strains})
    df = df.explode(column='P. aeruginosa Strain')

    columns = {'Regulator (TF or sigma)': 'regulator_locus_tag',
               'Target gene': 'gene_locus_tag',
               'mode of regulation': 'regulatory_effect',
               'Experimental Evidence': 'evidence',
               'PubMed Referencea': 'publication',
               'P. aeruginosa Strain': 'strain'}
    df = df.rename(columns=columns)
    df = df.assign(mechanism=None, effector_name=None)
    return df


class LiteratureTransformer(MultiStackTransformer, source='literature', version='0.0.0', register=False):
    _stack = ['faria_2016.xlsx',
              'fang_2017.xlsx',
              'turkarslan_2015.xls',
              'vasquez_2011.xls']
    _taxa = ['224308',
             '511145',
             '83332',
             None]
    _paer_strains = {'PAO1': '208964',
                     'PA103': '1081927',
                     'PA14': '652611',
                     'PAK': '1009714'}
    _source = ['bsub_faria_et_al_2017',
               'ecol_fang_et_al_2017',
               'mtub_turkarslan_et_al_2015',
               'paer_vasquez_et_al_2011']
    _net_reader = [read_bsub_faria_et_al_2017,
                   read_ecol_fang_et_al_2017,
                   read_mtub_turkarslan_et_al_2015,
                   read_paer_vasquez_et_al_2011]

    default_transform_stack = {
        'network': MultiStack(
            stack=_stack,
            taxa=_taxa,
            source=_source,
            reader=_net_reader
        ),
    }

    default_network_columns = SetList(['regulator_locus_tag', 'gene_locus_tag',
                                       'regulatory_effect', 'evidence', 'effector_name', 'mechanism',
                                       'publication', 'taxonomy', 'source'])

    def read_network(self) -> pd.DataFrame:
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)

        # special case of the paer_vasquez_et_al_2011 source which contains several strains in the same file,
        # thus having several organisms/taxonomy identifiers in the same file

        if 'strain' in network:
            taxonomy = network['strain'].map(self._paer_strains, na_action='ignore')
            taxonomy = network['taxonomy'].fillna(taxonomy)
            network = network.assign(taxonomy=taxonomy)

        filtered_networks = [self.filter_bsub_locus(network, 'regulator_locus_tag'),
                             self.filter_ecol_locus(network, 'regulator_locus_tag'),
                             self.filter_mtub_locus(network, 'regulator_locus_tag'),
                             self.filter_paer(network, 'regulator_locus_tag')]
        network = pd.concat(filtered_networks)
        network = network.reset_index(drop=True)
        return network

    def _transform_gene(self, network: pd.DataFrame, col: str) -> pd.DataFrame:
        network = network.assign(locus_tag=network[col].copy(), name=None)

        network = apply_processors(network, locus_tag=to_set_list)
        network = network.explode(column='locus_tag')

        network = apply_processors(network, locus_tag=[rstrip, lstrip])
        network = network.dropna(subset=['locus_tag'])
        network = drop_empty_string(network, 'locus_tag')
        network = drop_duplicates(df=network, subset=['locus_tag', 'taxonomy'], perfect_match=True)

        filtered_networks = [self.filter_bsub_locus(network, 'locus_tag'),
                             self.filter_ecol_locus(network, 'locus_tag'),
                             self.filter_mtub_locus(network, 'locus_tag'),
                             self.filter_paer_locus(network, 'locus_tag'),
                             self.filter_paer_names(network, 'locus_tag')]
        network = pd.concat(filtered_networks)
        network = network.reset_index(drop=True)

        gene_input_value = network['locus_tag'] + network['taxonomy']
        network = network.assign(input_value=gene_input_value)
        return network

    @staticmethod
    def merge_annotations(annotated: pd.DataFrame, original: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(annotated, original, on='input_value', suffixes=('_annotation', '_literature'))

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_literature')
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value'])
        return df

    @staticmethod
    def filter_bsub_locus(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df = df.copy()
        mask = (df[col].str.startswith('BSU')) & (df['source'] == 'bsub_faria_et_al_2017')
        df = df[mask]
        return df

    @staticmethod
    def filter_ecol_locus(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df = df.copy()
        mask = (df[col].str.startswith('b')) & (df['source'] == 'ecol_fang_et_al_2017')
        df = df[mask]
        return df

    @staticmethod
    def filter_mtub_locus(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df = df.copy()
        mask = (df[col].str.startswith('R')) & (df['source'] == 'mtub_turkarslan_et_al_2015')
        df = df[mask]
        return df

    @staticmethod
    def filter_paer(df: pd.DataFrame, _: str) -> pd.DataFrame:
        df = df.copy()
        mask = df['source'] == 'paer_vasquez_et_al_2011'
        df = df[mask]
        return df

    @staticmethod
    def filter_paer_locus(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df = df.copy()
        mask = (df[col].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df[mask]
        return df

    @staticmethod
    def filter_paer_names(df: pd.DataFrame, col: str) -> pd.DataFrame:
        # several loci in the paer_vasquez source are actually gene names,
        # so they should be reallocated to the name collum
        df = df.copy()
        mask = (~ df[col].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df[mask]
        df = df.rename(columns={col: 'name'})
        df = df.assign(**{col: None})
        return df

    @abstractmethod
    def transform(self):
        pass


class LiteratureConnector(Connector, source='literature', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
