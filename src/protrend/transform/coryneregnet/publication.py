import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.model import Publication, Regulator, Organism, Gene, TFBS, RegulatoryInteraction
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_str


class PublicationTransformer(CoryneRegNetTransformer,
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

    def transform_publication(self, network: pd.DataFrame) -> pd.DataFrame:
        pub = network.dropna(subset=['PMID'])
        pub = self.drop_empty_string(pub, 'PMID')
        pub = self.drop_duplicates(pub, subset=['PMID'])

        pub = pub.assign(pmid=pub['PMID'].copy())

        def split_pmid(item: str) -> list:
            items = item.split(',')
            return [item.lstrip().rstrip() for item in items]

        pub = apply_processors(pub, pmid=[to_str, split_pmid])

        pub = pub.explode('pmid')
        pub = pub.dropna(subset=['pmid'])
        pub = self.drop_empty_string(pub, 'pmid')
        pub = self.drop_duplicates(pub, subset=['pmid'])

        pub = self.create_input_value(pub, col='pmid')
        return pub

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)

        publications = self.transform_publication(network)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_coryneregnet'))
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
    default_connect_stack = {'publication': 'integrated_publication.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='organism',
                                    source_on='taxonomy', target_on='ncbi_taxonomy',
                                    source_processors={'taxonomy': [to_int_str]},
                                    target_processors={'ncbi_taxonomy': [to_int_str]})
        self.stack_json(df)


class PublicationToRegulatorConnector(CoryneRegNetConnector,
                                      source='coryneregnet',
                                      version='0.0.0',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):
    default_connect_stack = {'publication': 'integrated_publication.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='regulator',
                                    source_on='TF_locusTag', target_on='TF_locusTag')
        self.stack_json(df)


class PublicationToGeneConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):
    default_connect_stack = {'publication': 'integrated_publication.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='gene',
                                    source_on='TG_locusTag', target_on='TG_locusTag')
        self.stack_json(df)


class PublicationToTFBSConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):
    default_connect_stack = {'publication': 'integrated_publication.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='tfbs',
                                    source_on='PMID', target_on='PMID')
        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(CoryneRegNetConnector,
                                                  source='coryneregnet',
                                                  version='0.0.0',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='rin',
                                    source_on='PMID', target_on='PMID')
        self.stack_json(df)
