import pandas as pd

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Pathway, Regulator
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str, to_list_nan


class PathwayTransformer(RegPreciseTransformer,
                         source='regprecise',
                         version='0.0.0',
                         node=Pathway,
                         order=100,
                         register=True):
    default_transform_stack = {'pathway': 'Pathway.json'}
    columns = SetList(['protrend_id', 'name', 'kegg_pathways',
                       'pathway_id', 'url', 'regulog'])
    read_columns = SetList(['pathway_id', 'name', 'url', 'regulog'])

    def transform_pathway(self, pathway: pd.DataFrame) -> pd.DataFrame:
        # noinspection DuplicatedCode
        pathway = pathway.dropna(subset=['name'])
        pathway = self.drop_empty_string(pathway, 'name')
        pathway = self.drop_duplicates(df=pathway, subset=['name'])

        pathway = apply_processors(pathway, pathway_id=to_int_str, name=[rstrip, lstrip])

        # clean the unknown and NULL pathways
        mask = (pathway['name'] != 'unknown') & (pathway['name'] != 'NULL')
        pathway = pathway[mask]

        pathway = self.create_input_value(pathway, 'name')
        return pathway

    def transform(self):
        pathway = read_from_stack(stack=self.transform_stack, key='pathway',
                                  columns=self.read_columns, reader=read_json_lines)

        pathways = self.transform_pathway(pathway)
        annotated_pathways = self.annotate_pathways(pathways)

        df = self.merge_annotations_by_name(annotated_pathways, pathways)

        self.stack_transformed_nodes(df)
        return df


class PathwayToRegulatorConnector(RegPreciseConnector,
                                  source='regprecise',
                                  version='0.0.0',
                                  from_node=Pathway,
                                  to_node=Regulator,
                                  register=True):
    default_connect_stack = {'pathway': 'integrated_pathway.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='pathway',
                                                     target='regulator',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={'pathway_id': [to_int_str]},
                                                     target_processors={'pathway': [to_list_nan]})
        target_df = target_df.explode('pathway')
        target_df = apply_processors(target_df, pathway=to_int_str)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='pathway_id', target_on='pathway')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)
