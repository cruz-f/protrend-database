from typing import List

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import OrganismDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise.settings import OrganismSettings, OrganismToSource
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import DefaultTransformer


class OrganismTransformer(DefaultTransformer):
    default_settings = OrganismSettings
    columns = {'genome_id', 'name', 'taxonomy', 'url', 'regulon',
               'species', 'strain', 'family', 'phylum',
               'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
               'genbank_accession', 'genbank_ftp',
               'ncbi_assembly', 'assembly_accession'}

    def _read_genome(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('genome')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['genome_id', 'name', 'taxonomy', 'url', 'regulon'])

        return df

    def _transform_genome(self, genome: pd.DataFrame) -> pd.DataFrame:

        genome = self.drop_duplicates(df=genome, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=genome, col='name')

        genome['input_value'] = genome['name']

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

        genome = self._read_genome()
        genome = self._transform_genome(genome)

        names = list(genome['input_value'])

        organisms = self._transform_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        df['name'] = df['name_annotation']

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class OrganismToSourceConnector(DefaultConnector):
    default_settings = OrganismToSource

    def _read_organism(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('organism')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=OrganismTransformer.columns)

        return df

    def _read_source(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('source')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=SourceTransformer.columns)

        return df

    def connect(self):

        organism = self._read_organism()
        source = self._read_source()

        from_identifiers = organism['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=organism['url'].tolist(),
                      external_identifier=organism['genome_id'].tolist(),
                      key=['genome_id'] * size)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)
