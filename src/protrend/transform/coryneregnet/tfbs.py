import pandas as pd

from protrend.io import read_json_frame, read_from_stack, read_from_multi_stack
from protrend.model import TFBS
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.transform.coryneregnet.organism import OrganismTransformer
from protrend.utils import SetList, is_null, build_stack
from protrend.utils.processors import apply_processors


class TFBSTransformer(CoryneRegNetTransformer,
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

    def transform_tfbs(self, network: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # adding the organism protrend id
        tfbs = pd.merge(network, organism, on='taxonomy')

        tfbs = tfbs.assign(sequence=tfbs['Binding_site'].copy())

        # dropping nulls
        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = self.drop_empty_string(tfbs, 'sequence')

        # process multiple bs
        def binding_site_sequences(sequences: str) -> list:
            if is_null(sequences):
                return []
            return sequences.split(';')

        tfbs = apply_processors(tfbs, sequence=binding_site_sequences)
        tfbs = tfbs.explode(column='sequence')

        # drop duplicated interactions, namely regulator-target-binding site-organism
        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = self.drop_empty_string(tfbs, 'sequence')
        tfbs = self.drop_duplicates(df=tfbs,
                                    subset=['TF_locusTag', 'TG_locusTag', 'sequence', 'organism'],
                                    perfect_match=True)

        return tfbs

    @staticmethod
    def site_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.assign(strand=None, start=None, stop=None, length=tfbs['sequence'].str.len())
        return tfbs

    def transform_organism(self, organism_stack: dict) -> pd.DataFrame:
        organism = read_from_stack(stack=organism_stack,
                                   key='organism',
                                   columns=OrganismTransformer.columns,
                                   reader=read_json_frame)
        organism = self.select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'ncbi_taxonomy': 'taxonomy', 'protrend_id': 'organism'})
        return organism

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)

        organism_stack = build_stack(source=self.source, version=self.version,
                                     stack_to_load={'organism': 'integrated_organism.json'}, sa=False)
        organism = self.transform_organism(organism_stack)

        tfbs = self.transform_tfbs(network, organism)
        tfbs = self.site_coordinates(tfbs)
        tfbs = self.site_hash(tfbs)

        self.stack_transformed_nodes(tfbs)
        return tfbs
