import pandas as pd

from protrend.io import read_json_lines, read
from protrend.model import Organism
from protrend.transform.mix_ins import OrganismMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value, merge_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str


class OrganismTransformer(OrganismMixIn, RegPreciseTransformer,
                          source='regprecise',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession',
                       'genome_id', 'taxonomy', 'url', 'regulon'])

    @staticmethod
    def transform_organism(genome: pd.DataFrame) -> pd.DataFrame:
        genome = genome.dropna(subset=['name'])
        genome = drop_empty_string(genome, 'name')
        genome = drop_duplicates(df=genome, subset=['name'])

        genome = apply_processors(genome, genome_id=to_int_str, name=[rstrip, lstrip])

        genome = create_input_value(genome, col='name')
        return genome

    def transform(self):
        genome = read(source=self.source, version=self.version,
                      file='Genome.json', reader=read_json_lines,
                      default=pd.DataFrame(columns=['genome_id', 'name', 'taxonomy', 'url', 'regulon']))

        organisms = self.transform_organism(genome)
        annotated_organisms = self.annotate_organisms(organisms)

        df = pd.merge(annotated_organisms, organisms, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
