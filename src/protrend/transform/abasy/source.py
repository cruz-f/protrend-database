from protrend.model import Source, Organism, Regulator, Gene, RegulatoryInteraction
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList
from protrend.utils.constants import DATABASE


class SourceTransformer(SourceMixIn, AbasyTransformer,
                        source='abasy',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['abasy']
    type = [DATABASE]
    url = ['https://abasy.ccg.unam.mx']
    doi = ['10.1016/j.csbj.2020.05.015']
    authors = [['Juan M. Escorcia-Rodríguez, AndreasTauch, Julio A. Freyre-Gonzáleza']]
    description = ['Abasy Atlas v2.2: The most comprehensive and up-to-date inventory of meta-curated, historical '
                   'bacterial regulatory networks, their completeness and system-level characterization']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(AbasyConnector,
                                source='abasy',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self.create_connection(source='source', target='organism', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatorConnector(AbasyConnector,
                                 source='abasy',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self.create_connection(source='source', target='regulator', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToGeneConnector(AbasyConnector,
                            source='abasy',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self.create_connection(source='source', target='gene', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(AbasyConnector,
                                             source='abasy',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self.create_connection(source='source', target='rin', cardinality='one_to_many')
        self.stack_connections(df)
