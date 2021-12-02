import re
from collections import defaultdict
from statistics import mode, StatisticsError
from typing import List, Union

import numpy as np
import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import TFBS, Source, Organism
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList, is_null
from protrend.utils.processors import (apply_processors, remove_ellipsis, upper_case, to_list, flatten_set_list,
                                       to_int_str, to_str, to_set_list, operon_hash, site_hash)

regprecise_tfbs_pattern = re.compile(r'-\([0-9]+\)-')


class TFBSTransformer(RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=TFBS,
                      order=70,
                      register=True):
    default_transform_stack = {'tfbs': 'TFBS.json', 'gene': 'integrated_gene.json'}
    columns = SetList(['position', 'score', 'sequence', 'tfbs_id_old', 'tfbs_id', 'url',
                       'regulon', 'operon', 'gene', 'gene_protrend_id', 'gene_strand',
                       'gene_start', 'gene_old_locus_tag', 'length', 'strand', 'start', 'stop',
                       'site_hash', 'protrend_id'])
    read_columns = SetList(['position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon', 'operon', 'gene'])

    @staticmethod
    def _reduce_sequence(position, sequence):

        all_subgroups: List[Union[re.Match, None]] = list(re.finditer(regprecise_tfbs_pattern, sequence))
        all_subgroups += [None]

        sequences = []
        last_pos = position
        start = 0

        for subgroup in all_subgroups:

            if subgroup is None:
                seq = sequence[start:len(sequence)]
                pos = last_pos
                sequences.append((seq, pos))

            else:
                group = subgroup.group().replace('-(', '').replace(')-', '')
                end = subgroup.start()

                seq = sequence[start:end]
                pos = last_pos
                sequences.append((seq, pos))

                start = subgroup.end()

                length = int(group)
                last_pos = last_pos + len(seq) + length

        return sequences

    @staticmethod
    def _new_tfbs_identifier(tfbs_id: str, pos: int) -> str:

        _, *genes = tfbs_id.split('_')

        tfbs = [str(pos)] + list(genes)

        return '_'.join(tfbs)

    def _normalize_sequence(self, df: pd.DataFrame) -> pd.DataFrame:

        res = defaultdict(list)

        for _, row in df.iterrows():

            position = row['position']
            score = row['score']
            sequence = row['sequence']
            tfbs_id = row['tfbs_id']
            url = row['url']
            regulon = row['regulon']
            operon = row['operon']
            gene = row['gene']
            gene_id = row['gene_protrend_id']
            strand = row['gene_strand']
            start = row['gene_start']
            old_locus = row['gene_old_locus_tag']

            new_seqs = self._reduce_sequence(position=position, sequence=sequence)

            for seq, pos in new_seqs:
                new_tfbs_id = self._new_tfbs_identifier(tfbs_id, pos)

                res['position'].append(pos)
                res['score'].append(score)
                res['sequence'].append(seq)
                res['tfbs_id_old'].append(tfbs_id)
                res['tfbs_id'].append(new_tfbs_id)
                res['url'].append(url)
                res['regulon'].append(regulon)
                res['operon'].append(operon)
                res['gene'].append(gene)
                res['gene_protrend_id'].append(gene_id)
                res['gene_strand'].append(strand)
                res['gene_start'].append(start)
                res['gene_old_locus_tag'].append(old_locus)

        return pd.DataFrame(res)

    def _transform_tfbs(self, tfbs: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        # filter by sequence
        tfbs = tfbs.dropna(subset=['sequence'])

        # filter by gene
        tfbs = apply_processors(tfbs, gene=to_list)
        tfbs = tfbs.explode(column='gene')

        tfbs = pd.merge(tfbs, gene, left_on='gene', right_on='gene_old_locus_tag')

        aggr = {'gene': to_set_list,
                'gene_protrend_id': to_set_list,
                'gene_strand': to_set_list,
                'gene_start': to_set_list,
                'gene_old_locus_tag': to_set_list,
                'regulon': flatten_set_list,
                'operon': flatten_set_list}
        tfbs = self.group_by(df=tfbs, column='tfbs_id', aggregation=aggr)

        # filter by regulon, sequence and position
        tfbs = apply_processors(tfbs, regulon=to_list)
        tfbs = tfbs.explode(column='regulon')
        tfbs = self.drop_duplicates(df=tfbs, subset=['tfbs_id', 'sequence', 'regulon'],
                                    perfect_match=True, preserve_nan=False)

        tfbs = apply_processors(tfbs, sequence=[remove_ellipsis, upper_case])

        tfbs = self._normalize_sequence(tfbs)

        tfbs = self.drop_duplicates(df=tfbs, subset=['tfbs_id', 'sequence', 'regulon'],
                                    perfect_match=True, preserve_nan=False)

        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def strand_mode(item):

            if is_null(item):
                return None

            try:
                m = mode(item)

                if is_null(m):
                    return None

                return m

            except StatisticsError:

                for sub_item in item:
                    return sub_item

        def start_forward(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmin(x)

        def start_reverse(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmax(x)

        tfbs['length'] = tfbs['sequence'].str.len()

        tfbs['strand'] = tfbs['gene_strand'].map(strand_mode, na_action='ignore')
        forward = tfbs['strand'] == 'forward'
        reverse = tfbs['strand'] == 'reverse'

        tfbs['start'] = None
        tfbs['stop'] = None

        tfbs.loc[forward, 'start'] = tfbs.loc[forward, 'gene_start'].map(start_forward, na_action='ignore')
        tfbs.loc[forward, 'start'] = tfbs.loc[forward, 'start'] + tfbs.loc[forward, 'position']

        tfbs.loc[reverse, 'start'] = tfbs.loc[reverse, 'gene_start'].map(start_reverse, na_action='ignore')
        tfbs.loc[reverse, 'start'] = tfbs.loc[reverse, 'start'] - tfbs.loc[reverse, 'position']

        tfbs.loc[forward, 'stop'] = tfbs.loc[forward, 'start'] + tfbs.loc[forward, 'length']
        tfbs.loc[reverse, 'stop'] = tfbs.loc[reverse, 'start'] + tfbs.loc[reverse, 'length']

        strand_mask = (tfbs['strand'] != 'reverse') & (tfbs['strand'] != 'forward')
        tfbs.loc[strand_mask, 'strand'] = None

        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=self.read_columns, reader=read_json_lines)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'strand', 'start', 'locus_tag_old')
        gene = gene.rename(columns={'strand': 'gene_strand', 'start': 'gene_start',
                                    'locus_tag_old': 'gene_old_locus_tag', 'protrend_id': 'gene_protrend_id'})
        gene = gene.dropna(subset=['gene_old_locus_tag', 'gene_protrend_id'])
        gene = self.drop_duplicates(df=gene, subset=['gene_old_locus_tag', 'gene_protrend_id'], preserve_nan=False)

        df = self._transform_tfbs(tfbs=tfbs, gene=gene)
        df = self._tfbs_coordinates(df)

        # filter by site hash: length + strand + start + genes
        df2 = apply_processors(df,
                               sequence=[to_str, to_list],
                               length=[to_str, to_list],
                               strand=[to_str, to_list],
                               start=[to_str, to_list],
                               gene_protrend_id=[to_list, operon_hash, to_list])

        df['site_hash'] = df2['sequence'] + df2['length'] + df2['strand'] + df2['start'] + df2['gene_protrend_id']
        df = apply_processors(df, site_hash=site_hash)

        df = self.drop_duplicates(df=df, subset=['site_hash'], perfect_match=True)
        df = df.dropna(subset=['site_hash'])

        self.stack_transformed_nodes(df)

        return df


class TFBSToSourceConnector(RegPreciseConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=TFBS,
                            to_node=Source,
                            register=True):
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'source': 'integrated_source.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=tfbs['url'].tolist(),
                      external_identifier=tfbs['regulon'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class TFBSToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=TFBS,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = apply_processors(tfbs, regulon=to_int_str)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = apply_processors(regulator, regulon_id=to_int_str)

        merged = pd.merge(tfbs, regulator, left_on='regulon', right_on='regulon_id', suffixes=('_tfbs', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_tfbs', 'protrend_id_regulator', 'organism_protrend_id'])
        merged = merged.drop_duplicates(subset=['protrend_id_tfbs', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_tfbs'].tolist()
        to_identifiers = merged['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
