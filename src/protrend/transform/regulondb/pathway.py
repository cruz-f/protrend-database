from typing import List

import pandas as pd

from protrend.io import read_txt
from protrend.io.utils import read_from_stack
from protrend.model.model import Pathway
from protrend.transform.annotation import annotate_pathways
from protrend.transform.dto import PathwayDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.utils import SetList


class PathwayTransformer(RegulondbTransformer):
    default_node = Pathway
    default_transform_stack = {'groups': 'groups.txt'}
    default_order = 100
    columns = SetList(['synonyms', 'kegg_pathways', 'name', 'group_id', 'group_name', 'protrend_id'])
    read_columns = SetList(['group_id', 'group_name'])

    def _transform_pathway(self, pathway: pd.DataFrame) -> pd.DataFrame:
        pathway = pathway.dropna(subset=['group_name'])
        pathway = self.drop_duplicates(df=pathway, subset=['group_name'], perfect_match=True, preserve_nan=False)

        pathway = apply_processors(pathway, group_name=[rstrip, lstrip])

        self.create_input_value(pathway, 'group_name')
        return pathway

    @staticmethod
    def _transform_pathways(names: List[str]):
        dtos = [PathwayDTO(input_value=name) for name in names]
        annotate_pathways(dtos=dtos, names=names)

        # name: List[str]
        # synonyms: List[str]
        # kegg_pathways: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        pathway = read_from_stack(stack=self.transform_stack, file='groups', default_columns=self.read_columns,
                                  reader=read_txt, skiprows=29, names=self.read_columns)
        pathway = self._transform_pathway(pathway)

        names = pathway['input_value'].tolist()
        pathways = self._transform_pathways(names)

        df = pd.merge(pathways, pathway, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
