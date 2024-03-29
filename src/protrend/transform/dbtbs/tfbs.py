import pandas as pd

from protrend.io import read_json_lines
from protrend.io.utils import read_organism, read, read_gene, read_regulator
from protrend.model import TFBS
from protrend.report import ProtrendReporter
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dbtbs.organism import OrganismTransformer
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList, is_null
from protrend.utils.constants import REVERSE, FORWARD, UNKNOWN
from protrend.utils.processors import apply_processors, upper_case


class TFBSTransformer(TFBSMixIn, DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.4',
                      node=TFBS,
                      order=90,
                      register=True):
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'identifier', 'url', 'regulation', 'absolute_position', 'pubmed', 'tf', 'gene'])

    @staticmethod
    def transform_tfbs(tfbs: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.explode(column='url')
        tfbs = tfbs.explode(column='regulation')
        tfbs = tfbs.explode(column='absolute_position')
        tfbs = tfbs.explode(column='sequence')
        tfbs = tfbs.explode(column='tf')
        tfbs = tfbs.explode(column='gene')

        # filter duplicates, nan and empty strings
        tfbs = tfbs.dropna(subset=['identifier', 'sequence', 'absolute_position'])
        tfbs = drop_empty_string(tfbs, 'identifier', 'sequence', 'absolute_position')
        tfbs = drop_duplicates(df=tfbs, subset=['identifier', 'sequence', 'absolute_position'], perfect_match=True)

        # processing
        tfbs = apply_processors(tfbs, sequence=upper_case)

        # adding organism
        tfbs = tfbs.assign(organism=organism.loc[0, 'organism'])
        return tfbs

    @staticmethod
    def transform_organism(organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id')
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
                return UNKNOWN

            if isinstance(item, str):
                try:
                    seq_start, seq_stop = item.split('..')

                    seq_start = int(seq_start)
                    seq_stop = int(seq_stop)
                    diff = seq_stop - seq_start

                    if diff < 0:
                        return REVERSE

                    return FORWARD

                except ValueError:
                    return UNKNOWN

            return UNKNOWN

        length = tfbs['sequence'].map(sequence_len, na_action='ignore')
        start = tfbs['absolute_position'].map(sequence_start, na_action='ignore')
        stop = tfbs['absolute_position'].map(sequence_stop, na_action='ignore')
        strand = tfbs['absolute_position'].map(sequence_strand, na_action='ignore')

        tfbs = tfbs.assign(length=length, start=start, stop=stop, strand=strand)
        return tfbs

    def transform(self):
        tfbs = read(source=self.source, version=self.version, file='TFBS.json',
                    reader=read_json_lines,
                    default=pd.DataFrame(columns=['identifier', 'url', 'regulation', 'absolute_position', 'sequence',
                                                  'pubmed', 'tf', 'gene']))

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)

        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs, organism)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=tfbs.shape[0], properties=tfbs.shape[1])

        tfbs = self.site_coordinates(tfbs)

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        gene = gene.drop(columns=['tf'])

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        tfbs = pd.merge(tfbs, gene, left_on='gene', right_on='dbtbs_name')
        tfbs = pd.merge(tfbs, regulator, left_on='tf', right_on='dbtbs_name')

        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
