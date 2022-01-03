import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.model import Regulator
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip


class RegulatorTransformer(CoryneRegNetTransformer,
                           source='coryneregnet',
                           version='0.0.0',
                           node=Regulator,
                           order=100,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source'])

    def transform_regulator(self, network: pd.DataFrame) -> pd.DataFrame:
        regulator = network.dropna(subset=['TF_locusTag'])
        regulator = self.drop_empty_string(regulator, 'TF_locusTag')
        regulator = self.drop_duplicates(df=regulator, subset=['TF_locusTag'])

        regulator = apply_processors(regulator, TF_locusTag=[rstrip, lstrip], TF_name=[rstrip, lstrip])

        regulator = regulator.assign(locus_tag=regulator['TF_locusTag'].copy(), name=regulator['TF_name'].copy(),
                                     ncbi_taxonomy=regulator['taxonomy'].copy(), mechanism=None)

        mask = regulator['Is_sigma_factor'] == 'yes'
        regulator.loc[mask, 'mechanism'] = 'sigma factor'
        regulator.loc[~mask, 'mechanism'] = 'transcription factor'

        regulator = self.create_input_value(df=regulator, col='locus_tag')
        return regulator

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)

        regulators = self.transform_regulator(network)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

        self.stack_transformed_nodes(df)
        return df
