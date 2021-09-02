import re
from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.transform.connector import DefaultConnector
from protrend.transform.processors import (apply_processors, remove_ellipsis, upper_case, tfbs_left_position,
                                           operon_left_position, operon_strand, tfbs_right_position, null_to_none,
                                           to_int)
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.settings import TFBSSettings, TFBSToSource, TFBSToOrganism
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import Transformer

regprecise_tfbs_pattern = re.compile(r'-\([0-9]+\)-')


class TFBSTransformer(Transformer):
    default_settings = TFBSSettings
    columns = {'protrend_id',
               'position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon', 'operon', 'gene',
               'tfbs_id_old', 'position_left', 'position_right'}
    read_columns = {'position', 'score', 'sequence', 'tfbs_id', 'url', 'regulon', 'operon', 'gene'}

    @staticmethod
    def _reduce_sequence(position, sequence):

        all_subgroups: List[Union[re.Match, None]] = list(re.finditer(regprecise_tfbs_pattern, sequence))

        all_subgroups: List[Union[re.Match, None]] = all_subgroups + [None]

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

        return pd.DataFrame(res)

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:

        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = tfbs.explode('regulon')
        tfbs = self.drop_duplicates(df=tfbs, subset=['tfbs_id', 'sequence', 'regulon'],
                                    perfect_match=True, preserve_nan=False)
        tfbs = tfbs.reset_index(drop=True)

        apply_processors(remove_ellipsis, upper_case, df=tfbs, col='sequence')

        tfbs = self._normalize_sequence(tfbs)
        tfbs = self.drop_duplicates(df=tfbs, subset=['tfbs_id', 'sequence', 'regulon'],
                                    perfect_match=True, preserve_nan=False)
        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        gene = gene.set_index(gene['locus_tag_old'])
        apply_processors(null_to_none, df=gene, col='strand')
        apply_processors(null_to_none, df=gene, col='position_left')

        strands = []
        positions_left = []
        positions_right = []

        for position, seq, genes in zip(tfbs['position'], tfbs['sequence'], tfbs['gene']):

            op_strand = None

            for ge in genes:
                if ge in gene.index:
                    ge_strand = gene.loc[ge, 'strand']
                    op_strand = operon_strand(previous_strand=op_strand,
                                              current_strand=ge_strand)

            operon_left = None

            for ge in genes:
                if ge in gene.index:
                    ge_left = gene.loc[ge, 'position_left']
                    operon_left = operon_left_position(strand=op_strand,
                                                       previous_left=operon_left,
                                                       current_left=ge_left)

            tfbs_left = tfbs_left_position(strand=op_strand,
                                           gene_position=operon_left,
                                           gene_relative_position=position)

            tfbs_right = tfbs_right_position(strand=op_strand,
                                             gene_position=operon_left,
                                             gene_relative_position=position,
                                             tfbs_length=len(seq))

            strands.append(op_strand)
            positions_left.append(tfbs_left)
            positions_right.append(tfbs_right)

        tfbs['strand'] = strands
        tfbs['position_left'] = positions_left
        tfbs['position_right'] = positions_right

        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self._transform_stack, file='tfbs',
                               default_columns=self.read_columns, reader=read_json_lines)
        tfbs = self._transform_tfbs(tfbs)

        gene = read_from_stack(stack=self._transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'strand', 'position_left', 'locus_tag_old')
        gene = gene.dropna(subset=['locus_tag_old'])
        gene = self.drop_duplicates(df=gene, subset=['locus_tag_old'], perfect_match=False, preserve_nan=False)

        df = self._tfbs_coordinates(tfbs, gene)

        self._stack_transformed_nodes(df)

        return df


class TFBSToSourceConnector(DefaultConnector):
    default_settings = TFBSToSource

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

        self.stack_csv(df)


class TFBSToOrganismConnector(DefaultConnector):
    default_settings = TFBSToOrganism

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        merged = pd.merge(tfbs, regulator, left_on='regulon', right_on='regulon_id', suffixes=('_tfbs', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_tfbs'])
        merged = merged.dropna(subset=['protrend_id_regulator'])
        merged = merged.dropna(subset=['organism_protrend_id'])
        merged = merged.drop_duplicates(subset=['protrend_id_tfbs', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_tfbs'].tolist()
        to_identifiers = merged['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
