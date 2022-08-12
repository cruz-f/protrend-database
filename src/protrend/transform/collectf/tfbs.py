import pandas as pd

from protrend.io import read_json_lines
from protrend.io.utils import read_regulator, read
from protrend.model import TFBS
from protrend.report import ProtrendReporter
from protrend.transform.collectf.base import CollecTFTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.transformations import select_columns, drop_empty_string, group_by
from protrend.utils import SetList, is_null
from protrend.utils.constants import FORWARD, REVERSE, UNKNOWN
from protrend.utils.processors import (apply_processors, flatten_set_list_nan, to_set_list, to_list_nan, lstrip, rstrip,
                                       take_first)


class TFBSTransformer(TFBSMixIn, CollecTFTransformer,
                      source='collectf',
                      version='0.0.1',
                      node=TFBS,
                      order=80,
                      register=True):
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence',
                       'pubmed', 'regulon', 'operon', 'gene', 'experimental_evidence',
                       'regulon_id', 'regulator_name', 'tfbs', 'organism'])

    @staticmethod
    def transform_regulator(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'regulon_id', 'url', 'tfbs',
                                   'organism_protrend_id', 'organism_name', 'ncbi_taxonomy')

        regulator = apply_processors(regulator, tfbs=to_list_nan)
        regulator = regulator.explode('tfbs')
        regulator = regulator.dropna(subset=['tfbs'])
        regulator = drop_empty_string(regulator, 'tfbs')
        return regulator

    @staticmethod
    def transform_tfbs(tfbs: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.drop(columns=['organism'])

        tfbs = tfbs.rename(columns={'site_start': 'start', 'site_end': 'stop', 'site_strand': 'strand'})

        tfbs = apply_processors(tfbs, sequence=[rstrip, lstrip])
        # filter by sequence, start, stop, strand
        tfbs = tfbs.dropna(subset=['sequence', 'start', 'stop', 'strand'])
        tfbs = drop_empty_string(tfbs, 'sequence')

        tfbs = pd.merge(tfbs, regulator, left_on='tfbs_id', right_on='tfbs')

        tfbs = tfbs.drop_duplicates(subset=['tfbs_id', 'organism_protrend_id'])
        tfbs = drop_empty_string(tfbs, 'organism_protrend_id')
        tfbs = tfbs.dropna(subset=['organism_protrend_id'])

        tfbs = tfbs.rename(columns={'organism_protrend_id': 'organism'})
        return tfbs

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def strand_translation(item):
            if is_null(item):
                return

            if item == '1':
                return FORWARD

            elif item == '-1':
                return REVERSE

            return UNKNOWN

        def sequence_len(item):
            if is_null(item):
                return

            if isinstance(item, str):
                return len(item)

            return

        length = tfbs['sequence'].map(sequence_len, na_action='ignore')
        strand = tfbs['strand'].map(strand_translation, na_action='ignore')
        tfbs = tfbs.assign(length=length, strand=strand)
        return tfbs

    def transform(self):
        tfbs = read(source=self.source, version=self.version, file='TFBS.json',
                    reader=read_json_lines,
                    default=pd.DataFrame(
                        columns=['tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence',
                                 'pubmed', 'organism', 'regulon', 'operon', 'gene',
                                 'experimental_evidence'])
                    )

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)
        regulator = self.transform_regulator(regulator)

        tfbs = self.transform_tfbs(tfbs=tfbs, regulator=regulator)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=tfbs.shape[0], properties=tfbs.shape[1])

        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
