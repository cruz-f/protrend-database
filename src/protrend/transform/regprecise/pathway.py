from typing import List

import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.transform.annotation import annotate_pathways
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import PathwayDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, to_int_str
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.settings import PathwaySettings, PathwayToSource, PathwayToRegulator
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import Transformer


class PathwayTransformer(Transformer):
    default_settings = PathwaySettings
    columns = {'protrend_id',
               'pathway_id', 'name', 'url', 'regulog',
               'synonyms', 'kegg_pathways'}
    read_columns = {'pathway_id', 'name', 'url', 'regulog'}

    def _transform_pathway(self, pathway: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=pathway, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=df, col='name')
        apply_processors(to_int_str, df=df, col='pathway_id')
        apply_processors(to_int_str, df=df, col='regulog')

        self.create_input_value(df, 'name')
        return df

    @staticmethod
    def _transform_pathways(names: List[str]):
        dtos = [PathwayDTO(input_value=name) for name in names]
        annotate_pathways(dtos=dtos, names=names)

        # name: List[str]
        # synonyms: List[str]
        # kegg_pathways: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        pathway = read_from_stack(stack=self._transform_stack, file='pathway',
                                  default_columns=self.read_columns, reader=read_json_lines)
        pathway = self._transform_pathway(pathway)

        names = pathway['input_value'].tolist()
        pathways = self._transform_pathways(names)

        df = pd.merge(pathways, pathway, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df


class PathwayToSourceConnector(DefaultConnector):
    default_settings = PathwayToSource

    def connect(self):
        pathway = read_from_stack(stack=self._connect_stack, file='pathway',
                                  default_columns=PathwayTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = pathway['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=pathway['url'].tolist(),
                      external_identifier=pathway['pathway_id'].tolist(),
                      key=['pathway_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)


class PathwayToRegulatorConnector(DefaultConnector):
    default_settings = PathwayToRegulator

    def connect(self):
        pathway = read_from_stack(stack=self._connect_stack, file='pathway',
                                  default_columns=PathwayTransformer.columns, reader=read_json_frame)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.explode('pathway')

        merged = pd.merge(pathway, regulator, left_on='pathway_id', right_on='pathway',
                          suffixes=('_pathway', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_pathway'])
        merged = merged.dropna(subset=['protrend_id_regulator'])
        merged = merged.drop_duplicates(subset=['protrend_id_pathway', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_pathway'].tolist()
        to_identifiers = merged['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class PathwayToGeneConnector(DefaultConnector):
    default_settings = PathwayToRegulator

    def connect(self):
        pathway = read_from_stack(stack=self._connect_stack, file='pathway',
                                  default_columns=PathwayTransformer.columns, reader=read_json_frame)

        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        apply_processors(list, df=gene, col='regulon')
        gene = gene.explode('regulon')

        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        apply_processors(list, df=regulator, col='pathway')
        regulator = regulator.explode('pathway')

        merged = pd.merge(pathway, regulator, left_on='pathway_id', right_on='pathway',
                          suffixes=('_pathway', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_pathway'])
        merged = merged.dropna(subset=['protrend_id_regulator'])
        merged = merged.drop_duplicates(subset=['protrend_id_pathway', 'protrend_id_regulator'])

        merged = pd.merge(merged, gene, left_on='regulon_id', right_on='regulon',
                          suffixes=('_pathway_regulator', '_gene'))

        from_identifiers = merged['protrend_id_pathway'].tolist()
        to_identifiers = merged['protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
