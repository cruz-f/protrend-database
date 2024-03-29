import re
from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_organism, read_gene
from protrend.model import TFBS
from protrend.report import ProtrendReporter
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.transformations import select_columns, group_by, drop_empty_string
from protrend.utils import SetList
from protrend.utils.constants import FORWARD, REVERSE
from protrend.utils.processors import (apply_processors, remove_ellipsis, upper_case, flatten_set_list_nan,
                                       to_set_list, to_list_nan, take_last, strand_mode, start_forward, start_reverse,
                                       to_int_str)

regprecise_tfbs_pattern = re.compile(r'-\([0-9]+\)-')


class TFBSTransformer(TFBSMixIn, RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=TFBS,
                      order=70,
                      register=True):
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash'])

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

        if res:
            return pd.DataFrame(res)

        return pd.DataFrame(columns=['tfbs_id', 'position', 'sequence'])

    def transform_tfbs(self, tfbs: pd.DataFrame, gene: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # filter by sequence and position
        tfbs = tfbs.dropna(subset=['sequence', 'position'])
        tfbs = drop_empty_string(tfbs, 'sequence', 'position')

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
                'regulon': flatten_set_list_nan,
                'operon': flatten_set_list_nan,
                'gene': flatten_set_list_nan}
        tfbs = group_by(df=tfbs, column='tfbs_id', aggregation=aggr)

        ri = tfbs.drop(columns=['sequence', 'position'])

        sequences = select_columns(tfbs, 'tfbs_id', 'sequence', 'position')
        sequences = apply_processors(sequences, sequence=[remove_ellipsis, upper_case])
        sequences = self.transform_sequence(sequences)

        tfbs = pd.merge(ri, sequences, on='tfbs_id')

        # + 'organism', 'ncbi_taxonomy', 'genome'
        tfbs = pd.merge(tfbs, organism, on='ncbi_taxonomy')
        tfbs = tfbs.drop(columns=['genome'])
        return tfbs

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        strands = []
        starts = []
        stops = []
        lengths = []
        for strand, start, position, sequence in zip(tfbs['strand'], tfbs['start'], tfbs['position'], tfbs['sequence']):

            length = len(sequence)
            strand = strand_mode(strand)

            if strand == FORWARD:
                start = start_forward(start) + position
                stop = start + length

            elif strand == REVERSE:
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

    @staticmethod
    def transform_gene(gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'strand', 'start', 'tfbs', 'ncbi_taxonomy')
        gene = gene.rename(columns={'tfbs': 'tfbs_id'})
        gene = apply_processors(gene, tfbs_id=to_list_nan)
        gene = apply_processors(gene, ncbi_taxonomy=to_int_str)
        gene = gene.explode('tfbs_id')
        return gene

    def transform(self):
        tfbs = read(source=self.source, version=self.version,
                    file='TFBS.json', reader=read_json_lines,
                    default=pd.DataFrame(columns=['position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon',
                                                  'operon', 'gene']))
        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

        gene = self.transform_gene(gene)
        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs=tfbs, gene=gene, organism=organism)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=tfbs.shape[0], properties=tfbs.shape[1])

        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
