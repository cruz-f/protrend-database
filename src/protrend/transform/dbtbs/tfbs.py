import pandas as pd

from protrend.io import read_json_frame, read_json_lines, read_from_stack
from protrend.model import TFBS
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dbtbs.organism import OrganismTransformer
from protrend.utils import SetList, is_null
from protrend.utils.processors import apply_processors, upper_case


class TFBSTransformer(DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.4',
                      node=TFBS,
                      order=90,
                      register=True):
    default_transform_stack = {'tfbs': 'TFBS.json',
                               'organism': 'integrated_organism.json'}
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'identifier', 'url', 'regulation', 'absolute_position', 'pubmed', 'tf', 'gene'])
    read_columns = SetList(['identifier', 'url', 'regulation', 'absolute_position', 'sequence', 'pubmed', 'tf', 'gene'])

    def transform_tfbs(self, tfbs: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.explode(column='url')
        tfbs = tfbs.explode(column='regulation')
        tfbs = tfbs.explode(column='absolute_position')
        tfbs = tfbs.explode(column='sequence')
        tfbs = tfbs.explode(column='tf')
        tfbs = tfbs.explode(column='gene')

        # filter duplicates, nan and empty strings
        tfbs = tfbs.dropna(subset=['identifier', 'sequence', 'absolute_position'])
        tfbs = self.drop_empty_string(tfbs, 'identifier', 'sequence', 'absolute_position')
        tfbs = self.drop_duplicates(df=tfbs, subset=['identifier', 'sequence', 'absolute_position'], perfect_match=True)

        # processing
        tfbs = apply_processors(tfbs, sequence=upper_case)

        # adding organism
        tfbs = tfbs.reset_index(drop=True)
        organism = organism.reset_index(drop=True)
        tfbs = pd.concat([tfbs, organism], axis=1)
        return tfbs

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = self.select_columns(organism, 'protrend_id')
        organism = organism.rename(columns={'protrend_id': 'organism'})
        return organism

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def sequence_len(item):
            if is_null(item):
                return

            if isinstance(item, str):
                return len(item)

            return

        def sequence_start(item):
            if is_null(item):
                return

            if isinstance(item, str):
                try:
                    seq_start, _ = item.split('..')

                    seq_start = int(seq_start)
                    return seq_start

                except ValueError:
                    return

            return

        def sequence_stop(item):
            if is_null(item):
                return

            if isinstance(item, str):
                try:
                    _, seq_stop = item.split('..')

                    seq_stop = int(seq_stop)
                    return seq_stop

                except ValueError:
                    return

            return

        def sequence_strand(item):
            if is_null(item):
                return

            if isinstance(item, str):
                try:
                    seq_start, seq_stop = item.split('..')

                    seq_start = int(seq_start)
                    seq_stop = int(seq_stop)
                    diff = seq_stop - seq_start

                    if diff < 0:
                        return 'reverse'

                    return 'forward'

                except ValueError:
                    return

            return

        length = tfbs['sequence'].map(sequence_len, na_action='ignore')
        start = tfbs['absolute_position'].map(sequence_start, na_action='ignore')
        stop = tfbs['absolute_position'].map(sequence_stop, na_action='ignore')
        strand = tfbs['absolute_position'].map(sequence_strand, na_action='ignore')

        tfbs = tfbs.assign(length=length, start=start, stop=stop, strand=strand)
        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs', columns=self.read_columns,
                               reader=read_json_lines)

        # noinspection DuplicatedCode
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)

        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs, organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
