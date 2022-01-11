import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import TFBS
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.regulondb.base import RegulondbTransformer, regulondb_reader
from protrend.transform.regulondb.organism import OrganismTransformer
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors


class TFBSTransformer(TFBSMixIn, RegulondbTransformer,
                      source='regulondb',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):
    default_transform_stack = {'site': 'site.txt',
                               'organism': 'integrated_organism.json'}
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                       'site_internal_comment', 'key_id_org', 'site_length'])

    read_columns = SetList(['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                            'site_internal_comment', 'key_id_org', 'site_length'])

    @staticmethod
    def transform_tfbs(tfbs: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.assign(sequence=tfbs['site_sequence'].copy())

        # filter by nan and duplicates
        def remove_non_tfbs_sequence(sequence: str):

            tfbs_seq = ''
            for nucleotide in sequence:
                if nucleotide.isupper():
                    tfbs_seq += nucleotide

            if tfbs_seq:
                return tfbs_seq

            return

        tfbs = apply_processors(tfbs, sequence=remove_non_tfbs_sequence)

        tfbs = tfbs.dropna(subset=['site_id', 'sequence'])
        tfbs = drop_empty_string(tfbs, 'site_id', 'sequence')
        tfbs = drop_duplicates(df=tfbs, subset=['site_id', 'sequence'])

        # adding organism
        tfbs = tfbs.reset_index(drop=True)
        organism = organism.reset_index(drop=True)
        tfbs = pd.concat([tfbs, organism], axis=1)
        return tfbs

    @staticmethod
    def transform_organism(organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id')
        organism = organism.rename(columns={'protrend_id': 'organism'})
        return organism

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.assign(length=tfbs['sequence'].str.len(),
                           strand='forward',
                           start=tfbs['site_posleft'].copy(),
                           stop=tfbs['site_posright'].copy())
        return tfbs

    def transform(self):
        reader = regulondb_reader(skiprows=35, names=self.read_columns)
        tfbs = read_from_stack(stack=self.transform_stack, key='site', columns=self.read_columns,
                               reader=reader)

        # noinspection DuplicatedCode
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)

        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs, organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
