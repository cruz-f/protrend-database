import pandas as pd

from protrend.io.utils import read_organism, read
from protrend.model import TFBS
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, regulondb_reader
from protrend.transform.regulondb.organism import OrganismTransformer
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList
from protrend.utils.constants import FORWARD
from protrend.utils.processors import apply_processors


class TFBSTransformer(TFBSMixIn, RegulonDBTransformer,
                      source='regulondb',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
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
                           strand=FORWARD,
                           start=tfbs['site_posleft'].copy(),
                           stop=tfbs['site_posright'].copy())
        return tfbs

    def transform(self):
        columns = ['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                   'site_internal_comment', 'key_id_org', 'site_length']
        reader = regulondb_reader(skiprows=35, names=columns)
        tfbs = read(source=self.source, version=self.version,
                    file='site.txt', reader=reader,
                    default=pd.DataFrame(columns=columns))

        # noinspection DuplicatedCode
        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)

        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(tfbs, organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
