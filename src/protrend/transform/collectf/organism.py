import pandas as pd

from protrend.io import read_json_lines, read
from protrend.model import Organism, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.collectf.base import CollecTFTransformer, CollecTFConnector
from protrend.transform.mix_ins import OrganismMixIn
from protrend.transform.transformations import (merge_columns, create_input_value, drop_duplicates, group_by,
                                                drop_empty_string)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, take_last, flatten_set_list_nan


class OrganismTransformer(OrganismMixIn, CollecTFTransformer,
                          source='collectf',
                          version='0.0.1',
                          node=Organism,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession',
                       'genome_accession', 'taxonomy', 'regulon', 'tfbs'])

    @staticmethod
    def transform_organism(organism: pd.DataFrame) -> pd.DataFrame:
        organism = apply_processors(organism, name=[rstrip, lstrip])
        organism = organism.dropna(subset=['name'])
        organism = drop_empty_string(organism, 'name')

        aggregation = {'genome_accession': take_last, 'taxonomy': take_last}
        organism = group_by(df=organism, column='name', aggregation=aggregation, default=flatten_set_list_nan)

        organism = create_input_value(organism, col='name')
        return organism

    def transform(self):
        organism = read(source=self.source, version=self.version,
                        file='Organism.json', reader=read_json_lines,
                        default=pd.DataFrame(columns=['name', 'genome_accession', 'taxonomy', 'regulon', 'tfbs']))

        organisms = self.transform_organism(organism)
        annotated_organisms = self.annotate_organisms(organisms)

        df = pd.merge(annotated_organisms, organisms, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge name
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_collectf')

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class OrganismToRegulatorConnector(CollecTFConnector,
                                   source='collectf',
                                   version='0.0.1',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='regulator')
        self.stack_connections(df)


class OrganismToGeneConnector(CollecTFConnector,
                              source='collectf',
                              version='0.0.1',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='gene')
        self.stack_connections(df)


class OrganismToTFBSConnector(CollecTFConnector,
                              source='collectf',
                              version='0.0.1',
                              from_node=Organism,
                              to_node=TFBS,
                              register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='tfbs')
        self.stack_connections(df)


class OrganismToRegulatoryInteractionConnector(CollecTFConnector,
                                               source='collectf',
                                               version='0.0.1',
                                               from_node=Organism,
                                               to_node=RegulatoryInteraction,
                                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism')
        self.stack_connections(df)
