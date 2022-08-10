import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_organism
from protrend.model import Regulator
from protrend.transform.collectf.base import CollecTFTransformer
from protrend.transform.collectf.organism import OrganismTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import (drop_empty_string, drop_duplicates, create_input_value, select_columns,
                                                merge_columns, group_by)
from protrend.utils import SetList
from protrend.utils.constants import TRANSCRIPTION_FACTOR
from protrend.utils.processors import apply_processors, rstrip, lstrip, take_first, to_set_list, flatten_set_list_nan


class RegulatorTransformer(GeneMixIn, CollecTFTransformer,
                           source='collectf',
                           version='0.0.1',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop', 'mechanism',
                       'url', 'operon', 'gene', 'tfbs', 'experimental_evidence',
                       'organism_protrend_id', 'organism_name', 'ncbi_taxonomy',
                       'regulon_id', 'regulator_name'])

    def transform_regulon(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # ignore the uniprot annotation as it is incorrect
        regulon = regulon.drop(columns=['uniprot_accession'])

        # making a copy of the original name
        regulon = regulon.assign(regulator_name=regulon['name'].copy())

        regulon = apply_processors(regulon, name=[rstrip, lstrip], organism=[rstrip, lstrip])
        regulon = regulon.dropna(subset=['name', 'organism'])
        regulon = drop_empty_string(regulon, 'name', 'organism')

        df = pd.merge(regulon, organism, left_on='organism', right_on='organism_name')
        df = df.drop(columns=['organism'])

        # creating the regulon id: regulon_id = organism_name + '_' + regulator_name
        df['regulon_id'] = df.apply(lambda row: f'{row["organism_name"]}_{row["name"]}', axis=1)

        # group by regulon_id
        aggregation = {'name': take_first, 'url': to_set_list,
                       'operon': flatten_set_list_nan, 'gene': flatten_set_list_nan,
                       'tfbs': flatten_set_list_nan, 'experimental_evidence': flatten_set_list_nan,
                       'organism_protrend_id': take_first, 'organism_name': take_first, 'ncbi_taxonomy': take_first}
        df = group_by(df=df, column='regulon_id', aggregation=aggregation)

        df = df.assign(mechanism=TRANSCRIPTION_FACTOR)

        df = create_input_value(df=df, col='regulon_id')
        return df

    @staticmethod
    def transform_organism(organism: pd.DataFrame):
        organism = select_columns(organism, 'protrend_id', 'name', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism_protrend_id',
                                            'name': 'organism_name'})
        return organism

    def transform(self):
        regulon = read(source=self.source, version=self.version,
                       file='Regulon.json', reader=read_json_lines,
                       default=pd.DataFrame(columns=['uniprot_accession', 'name', 'url', 'organism', 'operon',
                                                     'gene', 'tfbs', 'experimental_evidence']))

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        organism = self.transform_organism(organism)

        regulators = self.transform_regulon(regulon, organism)
        annotated_regulators = self.annotate_genes(regulators)

        df = pd.merge(annotated_regulators, regulators, on='input_value', suffixes=('_annotation', '_collectf'))

        # drop regulators without a gene annotation, that is no locus tag
        df = df.dropna(subset=['locus_tag'])

        # merge name
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_collectf')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
