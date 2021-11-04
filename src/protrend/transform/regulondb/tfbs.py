import pandas as pd

from protrend.io import read_txt, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import TFBS
from protrend.transform.processors import (apply_processors, to_list, to_str, operon_hash, site_hash, to_set_list,
                                           take_first)
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class TFBSTransformer(RegulondbTransformer):
    default_node = TFBS
    default_transform_stack = {'site': 'site.txt',
                               'gene': 'integrated_gene.json'}
    default_order = 90
    columns = SetList(['sequence', 'strand', 'start', 'stop', 'length', 'site_hash', 'protrend_id',
                       'site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                       'site_internal_comment', 'key_id_org', 'site_length'])

    read_columns = SetList(['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                            'site_internal_comment', 'key_id_org', 'site_length'])

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        # filter by id
        tfbs = self.drop_duplicates(df=tfbs, subset=['site_id'], perfect_match=True, preserve_nan=True)

        # filter by nan
        tfbs = tfbs.dropna(subset=['site_id', 'site_sequence'])

        # filter out non-tfbs data
        def remove_non_tfbs_sequence(sequence: str):

            tfbs_seq = ''
            for nucleotide in sequence:
                if nucleotide.isupper():
                    tfbs_seq += nucleotide

            if tfbs_seq:
                return tfbs_seq

            return None

        tfbs = apply_processors(tfbs, site_sequence=remove_non_tfbs_sequence)

        # filter by coordinates
        tfbs = self.drop_duplicates(df=tfbs, subset=['site_sequence', 'site_posleft', 'site_posright'],
                                    perfect_match=True, preserve_nan=True)

        tfbs['sequence'] = tfbs['site_sequence']

        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def sequence_len(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                return len(item)

            return None

        tfbs['length'] = tfbs['sequence'].map(sequence_len, na_action='ignore')
        tfbs['start'] = tfbs['site_posleft']
        tfbs['stop'] = tfbs['site_posright']
        tfbs['strand'] = 'forward'

        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, file='site', default_columns=self.read_columns,
                               reader=read_txt, skiprows=35, names=self.read_columns)

        tfbs = self._transform_tfbs(tfbs)
        tfbs = self._tfbs_coordinates(tfbs)

        gene = read_from_stack(stack=self.transform_stack, file='gene', default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'gene_id', 'protrend_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        gene_tfbs = self._build_gene_tfbs()
        gene_tfbs = pd.merge(gene, gene_tfbs, on='gene_id')
        gene_tfbs = pd.merge(tfbs, gene_tfbs, on='site_id')
        aggregation = {'gene_id': to_set_list, 'gene_protrend_id': to_set_list}
        gene_tfbs = self.group_by(df=gene_tfbs, column='site_id', aggregation=aggregation, default=take_first)

        # filter by site hash: length + strand + start + genes
        df = apply_processors(gene_tfbs,
                              sequence=[to_str, to_list],
                              length=[to_str, to_list],
                              strand=[to_str, to_list],
                              start=[to_str, to_list],
                              gene_protrend_id=[to_list, operon_hash, to_list])

        gene_tfbs['site_hash'] = df['sequence'] + df['length'] + df['strand'] + df['start'] + df['gene_protrend_id']
        gene_tfbs = apply_processors(gene_tfbs, site_hash=site_hash)

        gene_tfbs = self.drop_duplicates(df=gene_tfbs, subset=['site_hash'], perfect_match=True, preserve_nan=True)
        gene_tfbs = gene_tfbs.dropna(subset=['site_hash'])

        self._stack_transformed_nodes(gene_tfbs)

        return gene_tfbs
