from typing import List

import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.annotation import annotate_organisms
from protrend.transform.connector import DefaultConnector
from protrend.transform.transformer import Transformer
from protrend.transform.dto import OrganismDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise.settings import OrganismSettings, OrganismToSource
from protrend.transform.regprecise import SourceTransformer


class OrganismTransformer(Transformer):
    default_settings = OrganismSettings
    columns = {'protrend_id',
               'genome_id', 'name', 'taxonomy', 'url', 'regulon',
               'species', 'strain', 'family', 'phylum',
               'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
               'genbank_accession', 'genbank_ftp',
               'ncbi_assembly', 'assembly_accession'}
    read_columns = {'genome_id', 'name', 'taxonomy', 'url', 'regulon'}

    def _transform_genome(self, genome: pd.DataFrame) -> pd.DataFrame:

        genome = self.drop_duplicates(df=genome, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=genome, col='name')

        genome = self.create_input_value(genome, col='name')

        return genome

    @staticmethod
    def _transform_organisms(names: List[str]):

        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, names=names)

        # name: List[str]
        # species: List[str]
        # strain: List[str]
        # family: List[str]
        # phylum: List[str]
        # ncbi_taxonomy: List[int]
        # refseq_accession: List[str]
        # refseq_ftp: List[str]
        # genbank_accession: List[str]
        # genbank_ftp: List[str]
        # ncbi_assembly: List[str]
        # assembly_accession: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):

        genome = read_from_stack(tl=self, file='genome', json=True, default_columns=self.read_columns)
        genome = self._transform_genome(genome)

        names = genome['input_value'].tolist()
        organisms = self._transform_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(['input_value'])

        self._stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(DefaultConnector):
    default_settings = OrganismToSource

    def connect(self):
        organism = read_from_stack(tl=self, file='organism', json=False, default_columns=OrganismTransformer.columns)
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

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

        self.stack_csv(df)
