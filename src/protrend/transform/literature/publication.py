import pandas as pd

from protrend.model import Publication, Organism, Regulator, Gene, RegulatoryInteraction
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector, read_literature_networks
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan, rstrip, lstrip


class PublicationTransformer(PublicationMixIn, LiteratureTransformer,
                             source='literature',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'evidence', 'effector_name', 'mechanism',
                       'publication', 'taxonomy', 'source'])

    @staticmethod
    def transform_publication(network: pd.DataFrame) -> pd.DataFrame:
        network = network.assign(pmid=network['publication'].copy())

        network = apply_processors(network, pmid=to_list_nan)
        network = network.explode(column='pmid')

        network = network.dropna(subset=['pmid'])
        network = drop_empty_string(network, 'pmid')
        network = drop_duplicates(df=network, subset=['pmid'], perfect_match=True)

        network = create_input_value(network, col='pmid')
        return network

    def transform(self):
        network = read_literature_networks(source=self.source, version=self.version)

        # noinspection DuplicatedCode
        publications = self.transform_publication(network)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_literature'))
        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(LiteratureConnector,
                                     source='literature',
                                     version='0.0.0',
                                     from_node=Publication,
                                     to_node=Organism,
                                     register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='organism',
                                    source_on='taxonomy', target_on='ncbi_taxonomy',
                                    source_processors={'taxonomy': [to_int_str]},
                                    target_processors={'ncbi_taxonomy': [to_int_str]})
        self.stack_connections(df)


class PublicationConnector(LiteratureConnector, register=False):

    def _connect(self, source: str, target: str) -> pd.DataFrame:
        source_df, target_df = self.transform_stacks(source=source,
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='publication',
                                                     target_on='publication',
                                                     source_processors={'publication': [to_list_nan]},
                                                     target_processors={'publication': [to_list_nan]})
        source_df = source_df.explode('publication')
        source_df = apply_processors(source_df, publication=[rstrip, lstrip])

        target_df = target_df.explode('publication')
        target_df = apply_processors(target_df, publication=[rstrip, lstrip])

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='publication', target_on='publication')

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids)


class PublicationToRegulatorConnector(PublicationConnector,
                                      source='literature',
                                      version='0.0.0',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):

    def connect(self):
        df = self._connect('publication', 'regulator')
        self.stack_connections(df)


class PublicationToGeneConnector(PublicationConnector,
                                 source='literature',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        df = self._connect('publication', 'gene')
        self.stack_connections(df)


class PublicationToRegulatoryInteractionConnector(PublicationConnector,
                                                  source='literature',
                                                  version='0.0.0',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):

    def connect(self):
        df = self._connect('publication', 'rin')
        self.stack_connections(df)
