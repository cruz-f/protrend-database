import numpy as np
import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import TFBS
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.transform.coryneregnet.gene import GeneTransformer
from protrend.utils.processors import (apply_processors, to_list, to_str, operon_hash, site_hash, to_set_list,
                                       take_first)
from protrend.utils import SetList, is_null


class TFBSTransformer(CoryneRegNetTransformer,
                      source='coryneregnet',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv',
                               'gene': 'integrated_gene.json'}
    columns = SetList(['sequence', 'strand', 'start', 'stop', 'length', 'site_hash', 'protrend_id',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role', 'TG_locusTag', 'TG_altLocusTag',
                       'TG_name', 'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy',
                       'Operon', 'Orientation', 'Genes', 'gene_protrend_id', 'gene_start'])

    def _transform_tfbs(self, regulation: pd.DataFrame) -> pd.DataFrame:
        # process multiple bs
        def binding_site_sequences(sequences: str) -> list:
            if is_null(sequences):
                return []
            return sequences.split(';')

        regulation = apply_processors(regulation, Binding_site=binding_site_sequences)
        regulation = regulation.explode(column='Binding_site')

        # filter by duplicated target gene - sequences
        regulation = self.drop_duplicates(df=regulation, subset=['TG_locusTag', 'Binding_site'], perfect_match=True)

        # filter by nan
        regulation = regulation.dropna(subset=['TG_locusTag', 'Binding_site'])

        regulation['sequence'] = regulation['Binding_site']

        return regulation

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def sequence_len(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                return len(item)

            return None

        def strand(item):
            if is_null(item):
                return None

            if item == '-':
                return 'reverse'

            if item == '+':
                return 'forward'

            return None

        def stop_forward(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmin(x)

        def stop_reverse(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmax(x)

        tfbs['length'] = tfbs['sequence'].map(sequence_len, na_action='ignore')

        tfbs['strand'] = tfbs['Orientation'].map(strand, na_action='ignore')
        forward = tfbs['strand'] == 'forward'
        reverse = tfbs['strand'] == 'reverse'

        tfbs['start'] = None
        tfbs['stop'] = None

        tfbs.loc[forward, 'stop'] = tfbs.loc[forward, 'gene_start'].map(stop_forward, na_action='ignore')
        tfbs.loc[reverse, 'stop'] = tfbs.loc[reverse, 'gene_start'].map(stop_reverse, na_action='ignore')

        tfbs.loc[forward, 'start'] = tfbs.loc[forward, 'stop'] - tfbs.loc[forward, 'length']
        tfbs.loc[reverse, 'start'] = tfbs.loc[reverse, 'stop'] + tfbs.loc[reverse, 'length']

        strand_mask = (tfbs['strand'] != 'reverse') & (tfbs['strand'] != 'forward')
        tfbs.loc[strand_mask, 'strand'] = None

        return tfbs

    def _transform_gene(self) -> pd.DataFrame:
        gene = read_from_stack(stack=self.transform_stack, file='gene', default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'start', 'stop', 'TG_locusTag')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'start': 'gene_start'})
        return gene

    def _transform_operon(self) -> pd.DataFrame:
        # 'Operon', 'Orientation', 'Genes'
        operon = self._build_operons()
        operon = operon.drop(columns='taxonomy')
        operon_by_gene = operon.explode(column='Genes')
        return operon_by_gene

    def _transform_operon_gene(self) -> pd.DataFrame:
        # 'Operon', 'Orientation', 'Genes'
        operon = self._transform_operon()
        # 'gene_protrend_id', 'TG_locusTag', 'gene_start'
        gene = self._transform_gene()
        operon_gene = pd.merge(operon, gene, left_on='Genes', right_on='TG_locusTag')

        aggregation = {'Orientation': take_first}
        operon_gene = self.group_by(df=operon_gene, column='Operon', aggregation=aggregation, default=to_set_list)
        return operon_gene

    def _transform_operon_tfbs(self) -> pd.DataFrame:
        # 'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role', 'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
        # 'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy', 'sequence'
        regulation = self._build_regulations()
        tfbs = self._transform_tfbs(regulation)
        tfbs = tfbs.drop(columns=['Operon'])

        # 'Operon', 'Orientation', 'Genes'
        operon = self._transform_operon()

        operon_tfbs = pd.merge(tfbs, operon, left_on='TG_locusTag', right_on='Genes')
        return operon_tfbs

    def transform(self):
        # 'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role', 'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
        # 'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy', 'sequence'
        # 'Operon', 'Orientation', 'Genes'
        operon_tfbs = self._transform_operon_tfbs()
        operon_tfbs = operon_tfbs.drop(columns=['Orientation', 'Genes'])

        # 'Operon', 'Orientation', 'Genes', 'gene_protrend_id', 'TG_locusTag', 'gene_start'
        operon_gene = self._transform_operon_gene()
        operon_gene = operon_gene.drop(columns=['TG_locusTag'])

        operon_gene_tfbs = pd.merge(operon_tfbs, operon_gene, on='Operon')
        operon_gene_tfbs = self._tfbs_coordinates(operon_gene_tfbs)

        # filter by site hash: length + strand + start + genes
        df = apply_processors(operon_gene_tfbs,
                              sequence=[to_str, to_list],
                              length=[to_str, to_list],
                              strand=[to_str, to_list],
                              start=[to_str, to_list],
                              gene_protrend_id=[to_list, operon_hash, to_list])

        operon_gene_tfbs['site_hash'] = df['sequence'] + df['length'] + df['strand'] + df['start'] + df[
            'gene_protrend_id']
        operon_gene_tfbs = apply_processors(operon_gene_tfbs, site_hash=site_hash)

        operon_gene_tfbs = self.drop_duplicates(df=operon_gene_tfbs, subset=['site_hash'], perfect_match=True)
        operon_gene_tfbs = operon_gene_tfbs.dropna(subset=['site_hash'])

        self.stack_transformed_nodes(operon_gene_tfbs)

        return operon_gene_tfbs
