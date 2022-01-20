import pandas as pd

from protrend.model import Publication, Regulator, Organism, Gene, TFBS, RegulatoryInteraction
from protrend.transform.coryneregnet.base import (CoryneRegNetTransformer, CoryneRegNetConnector,
                                                  read_coryneregnet_networks)
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value, merge_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_str


def _split_pmid(item: str) -> list:
    if ',' in item and ';' in item:
        sep = None

    elif ';' in item:
        sep = ';'

    elif ',' in item:
        sep = ','

    else:
        sep = None

    if sep is None:
        return []

    else:
        items = item.split(sep)
        return [item.lstrip().rstrip() for item in items]


class PublicationTransformer(PublicationMixIn, CoryneRegNetTransformer,
                             source='coryneregnet',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source'])

    @staticmethod
    def transform_publication(network: pd.DataFrame) -> pd.DataFrame:
        pub = network.dropna(subset=['PMID'])
        pub = drop_empty_string(pub, 'PMID')
        pub = drop_duplicates(pub, subset=['PMID'])

        pub = pub.assign(pmid=pub['PMID'].copy())

        pub = apply_processors(pub, pmid=[to_str, _split_pmid])

        pub = pub.explode('pmid')
        pub = pub.dropna(subset=['pmid'])
        pub = drop_empty_string(pub, 'pmid')
        pub = drop_duplicates(pub, subset=['pmid'])

        pub = create_input_value(pub, col='pmid')
        return pub

    def transform(self):
        network = read_coryneregnet_networks(self.source, self.version)

        publications = self.transform_publication(network)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_coryneregnet'))

        # merge pmid
        df = merge_columns(df=df, column='pmid', left='pmid_annotation', right='pmid_coryneregnet')

        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(CoryneRegNetConnector,
                                     source='coryneregnet',
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


class PublicationToRegulatorConnector(CoryneRegNetConnector,
                                      source='coryneregnet',
                                      version='0.0.0',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='regulator',
                                    source_on='TF_locusTag', target_on='TF_locusTag')
        self.stack_connections(df)


class PublicationToGeneConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='gene',
                                    source_on='TG_locusTag', target_on='TG_locusTag')
        self.stack_connections(df)


class PublicationToTFBSConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='tfbs',
                                    source_on='PMID', target_on='PMID')
        self.stack_connections(df)


class PublicationToRegulatoryInteractionConnector(CoryneRegNetConnector,
                                                  source='coryneregnet',
                                                  version='0.0.0',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='rin',
                                    source_on='PMID', target_on='PMID')
        self.stack_connections(df)
