import pandas as pd

from protrend.io import read_from_stack
from protrend.model import Publication, Regulator, TFBS, Gene, Organism, RegulatoryInteraction
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector, regulondb_reader
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList, build_stack
from protrend.utils.processors import apply_processors, to_int_str


class PublicationTransformer(PublicationMixIn, RegulonDBTransformer,
                             source='regulondb',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    default_transform_stack = {'publication': 'publication.txt'}
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'publication_id', 'reference_id', 'external_db_id', 'source'])
    read_columns = SetList(['publication_id', 'reference_id', 'external_db_id', 'author', 'title', 'source',
                            'year', 'publication_note', 'publication_internal_comment'])

    @staticmethod
    def transform_publication(publication: pd.DataFrame) -> pd.DataFrame:
        publication = select_columns(publication, 'publication_id', 'reference_id', 'external_db_id', 'source')
        publication = publication.assign(pmid=publication['reference_id'].copy())

        publication = apply_processors(df=publication, pmid=to_int_str)

        publication = publication.dropna(subset=['pmid'])
        publication = drop_empty_string(publication, 'pmid')
        publication = drop_duplicates(publication, subset=['pmid'])

        publication = create_input_value(publication, col='pmid')
        return publication

    def transform(self):
        reader = regulondb_reader(skiprows=36, names=self.read_columns)
        publication = read_from_stack(stack=self.transform_stack, key='publication', columns=self.read_columns,
                                      reader=reader)

        publications = self.transform_publication(publication)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = apply_processors(df, pmid=to_int_str)

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(RegulonDBConnector,
                                     source='regulondb',
                                     version='0.0.0',
                                     from_node=Publication,
                                     to_node=Organism,
                                     register=True):
    default_connect_stack = {'publication': 'integrated_publication.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='organism', cardinality='many_to_one')
        self.stack_json(df)


class PublicationConnector(RegulonDBConnector,
                           source='regulondb',
                           version='0.0.0',
                           register=False):
    obj_ev_pub = 'object_ev_method_pub_link.txt'

    def _connect(self, target: str, obj_ev_pub_col: str, target_col: str):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub_reader = regulondb_reader(skiprows=31, names=obj_ev_pub_cols)
        obj_ev_pub_stack = build_stack(self.source, self.version, stack_to_load={'obj_ev_pub': self.obj_ev_pub},
                                       dl=False)
        obj_ev_pub = read_from_stack(stack=obj_ev_pub_stack, key='obj_ev_pub',
                                     columns=obj_ev_pub_cols, reader=obj_ev_pub_reader)

        target_df = pd.merge(obj_ev_pub, target_df, left_on=obj_ev_pub_col, right_on=target_col)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='publication_id', target_on='publication_id')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        return df


class PublicationToRegulatorConnector(PublicationConnector,
                                      source='regulondb',
                                      version='0.0.0',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self._connect(target='regulator', obj_ev_pub_col='object_id', target_col='regulator_id')
        self.stack_json(df)


class PublicationToGeneConnector(PublicationConnector,
                                 source='regulondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'gene': 'integrated_gene.json'}

    def connect(self):
        df = self._connect(target='gene', obj_ev_pub_col='object_id', target_col='gene_id')
        self.stack_json(df)


class PublicationToTFBSConnector(PublicationConnector,
                                 source='regulondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='site_id')
        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(PublicationConnector,
                                                  source='regulondb',
                                                  version='0.0.0',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_tfbs.json'}

    def connect(self):
        regulator_df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='regulator_id')
        gene_df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='gene_id')
        tfbs_df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='site_id')
        df = pd.concat([regulator_df, gene_df, tfbs_df])
        df = df.reset_index(drop=True)
        self.stack_json(df)
