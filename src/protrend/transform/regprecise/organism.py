from typing import List

import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model.model import Organism, Source
from protrend.transform import OrganismDTO
from protrend.annotation import annotate_organisms
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class OrganismTransformer(RegPreciseTransformer):
    default_node = Organism
    default_transform_stack = {'genome': 'Genome.json'}
    default_order = 100
    columns = SetList(['species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly',
                       'assembly_accession', 'genome_id', 'taxonomy', 'url', 'regulon', 'name',
                       'protrend_id'])
    read_columns = SetList(['genome_id', 'name', 'taxonomy', 'url', 'regulon'])

    def _transform_genome(self, genome: pd.DataFrame) -> pd.DataFrame:
        genome = self.drop_duplicates(df=genome, subset=['name'], perfect_match=True, preserve_nan=False)

        genome = apply_processors(genome, name=[rstrip, lstrip], genome_id=to_int_str, taxonomy=to_int_str,
                                  regulon=to_int_str)

        genome = self.create_input_value(genome, col='name')

        return genome

    @staticmethod
    def _transform_organisms(names: List[str]):
        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, names=names)

        # name: List[str]
        # species: List[str]
        # strain: List[str]
        # ncbi_taxonomy: List[int]
        # refseq_accession: List[str]
        # refseq_ftp: List[str]
        # genbank_accession: List[str]
        # genbank_ftp: List[str]
        # ncbi_assembly: List[str]
        # assembly_accession: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        genome = read_from_stack(stack=self.transform_stack, file='genome',
                                 default_columns=self.read_columns, reader=read_json_lines)
        genome = self._transform_genome(genome)

        names = genome['input_value'].tolist()
        organisms = self._transform_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(RegPreciseConnector):
    default_from_node = Organism
    default_to_node = Source
    default_connect_stack = {'organism': 'integrated_organism.json', 'source': 'integrated_source.json'}

    def connect(self):
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = organism['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=organism['url'].tolist(),
                      external_identifier=organism['genome_id'].tolist(),
                      key=['genome_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
