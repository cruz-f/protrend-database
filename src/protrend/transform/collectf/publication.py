import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model import Publication, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan


class PublicationTransformer(CollectfTransformer,
                             source='collectf',
                             version='0.0.1',
                             node=Publication,
                             order=100,
                             register=True):
    default_transform_stack = {'tfbs': 'TFBS.json'}
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence',
                       'pubmed', 'organism', 'regulon', 'operon', 'gene', 'experimental_evidence'])
    read_columns = SetList(['tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence',
                            'pubmed', 'organism', 'regulon', 'operon', 'gene', 'experimental_evidence'])

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = apply_processors(df=tfbs, pubmed=to_list_nan)
        tfbs = tfbs.explode(column='pubmed')

        tfbs = tfbs.dropna(subset=['pubmed'])
        tfbs = self.drop_duplicates(df=tfbs, subset=['pubmed'])

        tfbs = tfbs.assign(pmid=tfbs['pubmed'].copy())

        tfbs = self.create_input_value(tfbs, col='pmid')
        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=self.read_columns, reader=read_json_lines)

        publications = self.transform_tfbs(tfbs)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_collectf'))
        df = apply_processors(df, pmid=to_int_str, year=to_int_str)
        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationConnector(CollectfConnector, register=False):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def _connect(self, target_column: str):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target='rin',
                                                     source_column='protrend_id',
                                                     target_column=target_column,
                                                     source_processors={'pmid': [to_int_str]},
                                                     target_processors={'pubmed': [to_list_nan]})
        target_df = target_df.explode('pubmed')
        target_df = apply_processors(target_df, pubmed=to_int_str)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='pmid', target_on='pubmed')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)


class PublicationToRegulatorConnector(PublicationConnector,
                                      source='collectf',
                                      version='0.0.1',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):

    def connect(self):
        return self._connect(target_column='regulator')


class PublicationToGeneConnector(PublicationConnector,
                                 source='collectf',
                                 version='0.0.1',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        return self._connect(target_column='gene')


class PublicationToTFBSConnector(PublicationConnector,
                                 source='collectf',
                                 version='0.0.1',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):

    def connect(self):
        return self._connect(target_column='tfbs')


class PublicationToRegulatoryInteractionConnector(PublicationConnector,
                                                  source='collectf',
                                                  version='0.0.1',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):

    def connect(self):
        return self._connect(target_column='protrend_id')
