import pandas as pd

from protrend.io import read_txt, read_json_frame, read_from_stack
from protrend.model.model import Promoter
from protrend.utils.processors import (apply_processors, to_list, to_str, operon_hash, upper_case, promoter_hash)
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.utils import SetList, is_null


class PromoterTransformer(RegulondbTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Promoter,
                          order=90,
                          register=True):
    default_transform_stack = {'promoter': 'promoter.txt',
                               'gene': 'integrated_gene.json'}
    columns = SetList(['sequence', 'strand', 'start', 'stop', 'length', 'promoter_hash', 'protrend_id',
                       'promoter_id', 'promoter_name', 'promoter_strand', 'pos_1', 'sigma_factor',
                       'basal_trans_val', 'equilibrium_const', 'kinetic_const', 'strength_seq',
                       'promoter_sequence', 'key_id_org', 'promoter_note', 'promoter_internal_comment'])

    read_columns = SetList(['promoter_id', 'promoter_name', 'promoter_strand', 'pos_1', 'sigma_factor',
                            'basal_trans_val', 'equilibrium_const', 'kinetic_const', 'strength_seq',
                            'promoter_sequence', 'key_id_org', 'promoter_note', 'promoter_internal_comment'])

    def _transform_promoter(self, promoter: pd.DataFrame) -> pd.DataFrame:
        # filter by id
        promoter = self.drop_duplicates(df=promoter, subset=['promoter_id'], perfect_match=True, preserve_nan=True)

        # filter by nan
        promoter = promoter.dropna(subset=['promoter_id', 'promoter_sequence'])

        promoter = apply_processors(promoter, promoter_sequence=upper_case)

        # filter by coordinates
        promoter = self.drop_duplicates(df=promoter,
                                        subset=['promoter_sequence', 'promoter_strand', 'pos_1'],
                                        perfect_match=True, preserve_nan=True)

        promoter['sequence'] = promoter['promoter_sequence']

        return promoter

    @staticmethod
    def _promoter_coordinates(promoter: pd.DataFrame) -> pd.DataFrame:

        def sequence_len(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                return len(item)

            return None

        def promoter_strand(item):
            if is_null(item):
                return None

            if isinstance(item, str):

                if item.lower() == 'reverse':
                    return 'reverse'

                elif item.lower() == 'forward':
                    return 'forward'

            return None

        promoter['length'] = promoter['sequence'].map(sequence_len, na_action='ignore')
        promoter['start'] = promoter['pos_1']
        promoter['stop'] = promoter['pos_1'] + promoter['length']
        promoter['strand'] = promoter['promoter_strand'].map(promoter_strand, na_action='ignore')

        return promoter

    def transform(self):
        promoter = read_from_stack(stack=self.transform_stack, file='promoter', default_columns=self.read_columns,
                                   reader=read_txt, skiprows=40, names=self.read_columns)

        promoter = self._transform_promoter(promoter)
        promoter = self._promoter_coordinates(promoter)

        gene = read_from_stack(stack=self.transform_stack, file='gene', default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'gene_id', 'protrend_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        gene_promoter = self._build_gene_promoter()
        gene_promoter = pd.merge(gene, gene_promoter, on='gene_id')

        promoter = pd.merge(promoter, gene_promoter, on='promoter_id')

        # filter by site hash: length + strand + start + genes
        df = apply_processors(promoter,
                              sequence=[to_str, to_list],
                              length=[to_str, to_list],
                              strand=[to_str, to_list],
                              start=[to_str, to_list],
                              gene_protrend_id=[to_list, operon_hash, to_list])

        _hash = df['sequence'] + df['length'] + df['strand'] + df['start'] + df['gene_protrend_id']
        promoter['promoter_hash'] = _hash
        promoter = apply_processors(promoter, promoter_hash=promoter_hash)

        promoter = self.drop_duplicates(df=promoter, subset=['promoter_hash'],
                                        perfect_match=True, preserve_nan=True)
        promoter = promoter.dropna(subset=['promoter_hash'])

        self._stack_transformed_nodes(promoter)

        return promoter
