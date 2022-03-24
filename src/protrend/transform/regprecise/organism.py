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
    manual_curation = {
        "Thermoanaerobacter ethanolicus X514":
            {"name": "Thermoanaerobacter sp. X514", "ncbi_taxonomy": 399726},
        "Marinobacter aqueolei":
            {"name": "Marinobacter nauticus", "ncbi_taxonomy": 2743},
        "Agrobacterium tumefaciens str. C58 (Cereon)":
            {"name": "Agrobacterium fabrum str. C58", "ncbi_taxonomy": 176299},
        "Silicibacter TM1040":
            {"name": "Ruegeria sp. TM1040", "ncbi_taxonomy": 292414},
        "Streptomyces coelicolor A3(2)":
            {"name": "Streptomyces coelicolor A3(2)", "ncbi_taxonomy": 100226},
        "Desulfovibrio desulfuricans G20":
            {"name": "Oleidesulfovibrio alaskensis G20", "ncbi_taxonomy": 207559},
        "Geobacter uraniumreducens Rf4":
            {"name": "Geotalea uraniireducens Rf4", "ncbi_taxonomy": 351605},
        "Blattabacterium sp. (Blattella germanica) str. Bge":
            {"name": "Blattabacterium sp. (Blattella germanica) str. Bge", "ncbi_taxonomy": 331104},
        "Lactobacillus rhamnosus GG":
            {"name": "Lacticaseibacillus rhamnosus GG", "ncbi_taxonomy": 568703},
        "Alteromonas macleodii 'Deep ecotype'":
            {"name": "Alteromonas mediterranea DE", "ncbi_taxonomy": 1774373}
    }

    def transform_organism(self, genome: pd.DataFrame) -> pd.DataFrame:
        genome = genome.dropna(subset=['name'])
        genome = drop_empty_string(genome, 'name')
        genome = drop_duplicates(df=genome, subset=['name'])

        genome = apply_processors(genome, genome_id=to_int_str, name=[rstrip, lstrip])

        def _map_name(key):
            map_dict = self.manual_curation.get(key, {'name': key})
            return map_dict['name']

        def _map_ncbi(key):
            map_dict = self.manual_curation.get(key, {'ncbi_taxonomy': None})
            return map_dict['ncbi_taxonomy']

        curated_names = genome['name'].map(_map_name)
        curated_tax = genome['name'].map(_map_ncbi)

        genome = genome.assign(name=curated_names, ncbi_taxonomy=curated_tax)
        genome = apply_processors(genome, name=[rstrip, lstrip], ncbi_taxonomy=to_int_str)

        genome = create_input_value(genome, col='name')
        return genome

    def transform(self):
        genome = read(source=self.source, version=self.version,
                      file='Genome.json', reader=read_json_lines,
                      default=pd.DataFrame(columns=['genome_id', 'name', 'taxonomy', 'url', 'regulon']))

        organisms = self.transform_organism(genome)
        annotated_organisms = self.annotate_organisms(organisms)

        df = pd.merge(annotated_organisms, organisms, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='ncbi_taxonomy', left='ncbi_taxonomy_annotation',
                           right='ncbi_taxonomy_regprecise')
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
