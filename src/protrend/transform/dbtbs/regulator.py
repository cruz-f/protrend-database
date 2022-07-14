import pandas as pd

from protrend.io import read_json_lines, read, read_genbank
from protrend.model import Regulator
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.constants import UNKNOWN, SIGMA_FACTOR, TRANSCRIPTION_FACTOR
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_list_nan, take_first


class RegulatorTransformer(GeneMixIn, DBTBSTransformer,
                           source='dbtbs',
                           version='0.0.4',
                           node=Regulator,
                           order=100,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop', 'mechanism',
                       'family', 'domain', 'domain_description', 'description', 'url',
                       'type', 'consensus_sequence', 'comment', 'subti_list', 'gene', 'tfbs',
                       'dbtbs_name'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        tf = tf.assign(dbtbs_name=tf['name'].copy())

        tf = apply_processors(tf, name=[rstrip, lstrip], family=[to_list_nan, take_first, rstrip, lstrip])

        # filter nan and duplicates
        tf = tf.dropna(subset=['name'])
        tf = drop_empty_string(tf, 'name')
        tf = drop_duplicates(df=tf, subset=['name'])

        tf = tf.assign(name_lower=tf['name'].str.lower(), mechanism=UNKNOWN)

        # sigma factors
        tf_sigma_mask = tf['family'] == 'Sigma factors'
        tf.loc[tf_sigma_mask, 'mechanism'] = SIGMA_FACTOR

        # tf mechanism
        tf_not_sigma_mask = tf['family'] != 'Sigma factors'
        tf.loc[tf_not_sigma_mask, 'mechanism'] = TRANSCRIPTION_FACTOR

        tf = pd.merge(tf, sequence, on='name_lower')
        tf = tf.drop(columns=['name_lower'])

        # for locus tag annotation
        tf = tf.assign(taxonomy='224308')

        tf = create_input_value(df=tf, col='locus_tag')
        return tf

    def transform(self):
        tf = read(source=self.source, version=self.version, file='TranscriptionFactor.json', reader=read_json_lines,
                  default=pd.DataFrame(columns=['name', 'family', 'domain', 'domain_description', 'description', 'url',
                                                'type', 'consensus_sequence', 'comment', 'subti_list', 'gene', 'tfbs']))

        sequence = read(source=self.source, version=self.version, file='protein_sequence.gb', reader=read_genbank,
                        default=pd.DataFrame(columns=['name_lower', 'locus_tag', 'genbank_accession',
                                                      'uniprot_accession']))

        regulators = self.transform_tf(tf=tf, sequence=sequence)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
