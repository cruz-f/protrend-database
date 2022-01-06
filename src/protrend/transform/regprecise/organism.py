import pandas as pd

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Organism
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str


class OrganismTransformer(RegPreciseTransformer,
                          source='regprecise',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
    default_transform_stack = {'genome': 'Genome.json'}
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession',
                       'genome_id', 'taxonomy', 'url', 'regulon'])
    read_columns = SetList(['genome_id', 'name', 'taxonomy', 'url', 'regulon'])

    def transform_organism(self, genome: pd.DataFrame) -> pd.DataFrame:
        genome = genome.dropna(subset=['name'])
        genome = self.drop_empty_string(genome, 'name')
        genome = self.drop_duplicates(df=genome, subset=['name'])

        genome = apply_processors(genome, genome_id=to_int_str, name=[rstrip, lstrip])

        genome = self.create_input_value(genome, col='name')
        return genome

    def transform(self):
        genome = read_from_stack(stack=self.transform_stack, key='genome',
                                 columns=self.read_columns, reader=read_json_lines)

        organisms = self.transform_organism(genome)
        annotated_organisms = self.annotate_organisms(organisms)

        df = pd.merge(annotated_organisms, organisms, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
