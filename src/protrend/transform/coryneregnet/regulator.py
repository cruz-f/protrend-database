import pandas as pd

from protrend.model import Regulator
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, read_coryneregnet_networks
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.constants import TRANSCRIPTION_FACTOR, SIGMA_FACTOR
from protrend.utils.processors import apply_processors, rstrip, lstrip


class RegulatorTransformer(GeneMixIn, CoryneRegNetTransformer,
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

    @staticmethod
    def transform_regulator(network: pd.DataFrame) -> pd.DataFrame:
        regulator = network.dropna(subset=['TF_locusTag'])
        regulator = drop_empty_string(regulator, 'TF_locusTag')
        regulator = drop_duplicates(df=regulator, subset=['TF_locusTag'])

        regulator = apply_processors(regulator, TF_locusTag=[rstrip, lstrip], TF_name=[rstrip, lstrip])

        regulator = regulator.assign(locus_tag=regulator['TF_locusTag'].copy(), name=regulator['TF_name'].copy(),
                                     ncbi_taxonomy=regulator['taxonomy'].copy(), mechanism=None)

        mask = regulator['Is_sigma_factor'] == 'yes'
        regulator.loc[mask, 'mechanism'] = SIGMA_FACTOR
        regulator.loc[~mask, 'mechanism'] = TRANSCRIPTION_FACTOR

        regulator = create_input_value(df=regulator, col='locus_tag')
        return regulator

    def transform(self):
        network = read_coryneregnet_networks(self.source, self.version)

        regulators = self.transform_regulator(network)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

        self.stack_transformed_nodes(df)
        return df
