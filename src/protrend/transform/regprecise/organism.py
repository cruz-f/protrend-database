from typing import List

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source
from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.dto import OrganismDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise.settings import OrganismSettings
from protrend.transform.transformer import Transformer


class OrganismTransformer(Transformer):

    def __init__(self, settings: OrganismSettings = None):

        if not settings:
            settings = OrganismSettings()

        super().__init__(settings)

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

        organisms = pd.DataFrame([dto.to_dict() for dto in dtos])

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

        if organisms.empty:
            organisms = pd.DataFrame(columns=['input_value', 'name',
                                              'species', 'strain', 'family', 'phylum',
                                              'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                                              'genbank_accession', 'genbank_ftp',
                                              'ncbi_assembly', 'assembly_accession'])

        return organisms

    def transform(self):

        genome = self._read_genome()
        genome = self._transform_genome(genome)

        names = list(genome['input_value'])

        organisms = self._transform_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        df['name'] = df['name_annotation']

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        from_identifiers = from_df['protrend_id'].tolist()

        size = len(from_identifiers)

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]
        to_identifiers = [regprecise_id] * size

        kwargs = dict(url=from_df['url'].tolist(),
                      external_identifier=from_df['genome_id'].tolist(),
                      key=['genome_id'] * size)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    kwargs=kwargs)

    def connect(self):

        source_connection = self._connect_to_source()
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, source_connection)
