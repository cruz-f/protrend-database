import pandas as pd

from protrend.io import read_json_lines
from protrend.io.utils import read_regulator, read
from protrend.model import TFBS
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
                       'pubmed', 'organism', 'regulon', 'operon', 'gene', 'experimental_evidence',
                       'uniprot_accession'])

    @staticmethod
    def transform_regulator(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'uniprot_accession', 'organism_protrend_id')
        return regulator

    @staticmethod
    def transform_tfbs(tfbs: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.rename(columns={'site_start': 'start', 'site_end': 'stop', 'site_strand': 'strand',
                                    'organism': 'organism_collectf'})

        tfbs = apply_processors(tfbs, sequence=[rstrip, lstrip], regulon=to_list_nan)
        # filter by sequence, start, stop, strand
        tfbs = tfbs.dropna(subset=['sequence', 'start', 'stop', 'strand'])
        tfbs = drop_empty_string(tfbs, 'sequence')

        # filter by organism
        tfbs = tfbs.explode(column='regulon')

        tfbs = pd.merge(tfbs, regulator, left_on='regulon', right_on='uniprot_accession')

        aggr = {'pubmed': flatten_set_list_nan, 'regulon': to_set_list, 'operon': flatten_set_list_nan,
                'experimental_evidence': flatten_set_list_nan, 'gene': flatten_set_list_nan}
        tfbs = group_by(df=tfbs, column='tfbs_id', aggregation=aggr, default=take_first)
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
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
