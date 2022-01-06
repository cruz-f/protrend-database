import re
from collections import defaultdict
from statistics import mode, StatisticsError
from typing import List, Union

import numpy as np
import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import TFBS
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.utils import SetList, is_null
from protrend.utils.processors import (apply_processors, remove_ellipsis, upper_case, to_list, flatten_set_list,
                                       to_set_list, to_list_nan, take_last)

regprecise_tfbs_pattern = re.compile(r'-\([0-9]+\)-')


class TFBSTransformer(RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=TFBS,
                      order=70,
                      register=True):
    default_transform_stack = {'tfbs': 'TFBS.json',
                               'organism': 'integrated_organism.json', 'gene': 'integrated_gene.json'}
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       ])
    read_columns = SetList(['position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon', 'operon', 'gene'])

    @staticmethod
    def split_sequence(position, sequence):

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

    def transform_sequence(self, df: pd.DataFrame) -> pd.DataFrame:

        res = defaultdict(list)

        for tfbs_id, position, sequence in zip(df['tfbs_id'], df['position'], df['sequence']):

            new_sequences = self.split_sequence(position=position, sequence=sequence)

            for new_sequence, new_position in new_sequences:
                res['tfbs_id'].append(tfbs_id)
                res['position'].append(new_position)
                res['sequence'].append(new_sequence)

        return pd.DataFrame(res)

    def transform_tfbs(self, tfbs: pd.DataFrame, gene: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # filter by sequence and position
        tfbs = tfbs.dropna(subset=['sequence', 'position'])
        tfbs = self.drop_empty_string(tfbs, 'sequence', 'position')

        # + 'strand', 'start', 'ncbi_taxonomy'
        tfbs = pd.merge(tfbs, gene, on='tfbs_id')

        # 'position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon', 'operon', 'gene'
        aggr = {'strand': to_set_list,
                'start': to_set_list,
                'ncbi_taxonomy': take_last,
                'position': take_last,
                'score': take_last,
                'sequence': take_last,
                'url': take_last,
                'regulon': flatten_set_list,
                'operon': flatten_set_list,
                'gene': flatten_set_list}
        tfbs = self.group_by(df=tfbs, column='tfbs_id', aggregation=aggr)

        ris = tfbs.drop(columns=['sequence', 'position'])

        sequences = self.select_columns(tfbs, 'tfbs_id', 'sequence', 'position')
        sequences = apply_processors(sequences, sequence=[remove_ellipsis, upper_case])
        sequences = self.transform_sequence(sequences)

        tfbs = pd.merge(ris, sequences, on='tfbs_id')

        # + 'organism', 'ncbi_taxonomy', 'genome'
        tfbs = pd.merge(tfbs, organism, on='ncbi_taxonomy')
        tfbs = tfbs.drop(columns=['genome'])
        return tfbs

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def strand_mode(item):

            if is_null(item):
                return

            try:
                m = mode(item)

                if is_null(m):
                    return

                return m

            except StatisticsError:

                for sub_item in item:
                    return sub_item

        def start_forward(item):
            if is_null(item):
                return

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmin(x)

        def start_reverse(item):
            if is_null(item):
                return

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmax(x)

        strands = []
        starts = []
        stops = []
        lengths = []
        for strand, start, position, sequence in zip(tfbs['strand'], tfbs['start'], tfbs['position'], tfbs['sequence']):

            length = len(sequence)
            strand = strand_mode(strand)

            if strand == 'forward':
                start = start_forward(start) + position
                stop = start + length

            elif strand == 'reverse':
                start = start_reverse(start) - position
                stop = start + length

            else:
                start = None
                stop = None

            strands.append(strand)
            starts.append(start)
            stops.append(stop)
            lengths.append(length)

        tfbs = tfbs.assign(strand=strands, start=starts, stop=stops, length=lengths)
        return tfbs

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'strand', 'start', 'tfbs', 'ncbi_taxonomy')
        gene = gene.rename(columns={'tfbs': 'tfbs_id'})
        gene = apply_processors(gene, tfbs_id=to_list_nan)
        gene = gene.explode('tfbs_id')
        return gene

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=self.read_columns, reader=read_json_lines)
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.columns, reader=read_json_frame)

        gene = self.transform_gene(gene)
        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs=tfbs, gene=gene, organism=organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
