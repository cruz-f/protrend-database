import pandas as pd

from protrend.io import read_json_lines, read
from protrend.model import Pathway, Regulator
from protrend.transform.mix_ins import PathwayMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str, to_list_nan


class PathwayTransformer(PathwayMixIn, RegPreciseTransformer,
                         source='regprecise',
                         version='0.0.0',
                         node=Pathway,
                         order=100,
                         register=True):
    columns = SetList(['protrend_id', 'name', 'kegg_pathways',
                       'pathway_id', 'url', 'regulog'])

    @staticmethod
    def transform_pathway(pathway: pd.DataFrame) -> pd.DataFrame:
        # noinspection DuplicatedCode
        pathway = pathway.dropna(subset=['name'])
        pathway = drop_empty_string(pathway, 'name')
        pathway = drop_duplicates(df=pathway, subset=['name'])

        pathway = apply_processors(pathway, pathway_id=to_int_str, name=[rstrip, lstrip])

        # clean the unknown and NULL pathways
        mask = (pathway['name'] != 'unknown') & (pathway['name'] != 'NULL')
        pathway = pathway[mask]

        pathway = create_input_value(pathway, 'name')
        return pathway

    def transform(self):
        pathway = read(source=self.source, version=self.version,
                       file='Pathway.json', reader=read_json_lines,
                       default=pd.DataFrame(columns=['pathway_id', 'name', 'url', 'regulog']))

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

    def connect(self):
        source_df, target_df = self.transform_stacks(source='pathway',
                                                     target='regulator',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='pathway_id',
                                                     target_on='pathway',
                                                     source_processors={'pathway_id': [to_int_str]},
                                                     target_processors={'pathway': [to_list_nan]})
        target_df = target_df.explode('pathway')
        target_df = apply_processors(target_df, pathway=to_int_str)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='pathway_id', target_on='pathway')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
