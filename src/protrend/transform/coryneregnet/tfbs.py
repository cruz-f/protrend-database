import pandas as pd

from protrend.io.utils import read_organism
from protrend.model import TFBS
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, read_coryneregnet_networks
from protrend.transform.coryneregnet.organism import OrganismTransformer
from protrend.transform.mix_ins import TFBSMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList, is_null
from protrend.utils.constants import UNKNOWN
from protrend.utils.processors import apply_processors, to_int_str


class TFBSTransformer(TFBSMixIn, CoryneRegNetTransformer,
                      source='coryneregnet',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):
    columns = SetList(['protrend_id', 'organism', 'start', 'stop', 'strand', 'sequence', 'length', 'site_hash',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source', 'organism'])

    @staticmethod
    def transform_tfbs(network: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, taxonomy=to_int_str)

        # adding the organism protrend id
        tfbs = pd.merge(network, organism, on='taxonomy')

        tfbs = tfbs.assign(sequence=tfbs['Binding_site'].copy())

        # dropping nulls
        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = drop_empty_string(tfbs, 'sequence')

        # process multiple bs
        def binding_site_sequences(sequences: str) -> list:
            if is_null(sequences):
                return []
            return sequences.split(';')

        tfbs = apply_processors(tfbs, sequence=binding_site_sequences)
        tfbs = tfbs.explode(column='sequence')

        # drop duplicated interactions, namely regulator-target-binding site-organism
        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = drop_empty_string(tfbs, 'sequence')
        tfbs = drop_duplicates(df=tfbs,
                               subset=['TF_locusTag', 'TG_locusTag', 'sequence', 'organism'],
                               perfect_match=True)

        return tfbs

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.assign(strand=UNKNOWN, start=None, stop=None, length=tfbs['sequence'].str.len())
        return tfbs

    @staticmethod
    def transform_organism(organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'ncbi_taxonomy': 'taxonomy', 'protrend_id': 'organism'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform(self):
        network = read_coryneregnet_networks(self.source, self.version)

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        organism = self.transform_organism(organism)

        tfbs = self.transform_tfbs(network, organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
