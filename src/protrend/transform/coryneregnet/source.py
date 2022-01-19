from protrend.model import Source, Organism, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList
from protrend.utils.constants import DATABASE


class SourceTransformer(SourceMixIn, CoryneRegNetTransformer,
                        source='coryneregnet',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['coryneregnet']
    type = [DATABASE]
    url = ['https://www.exbio.wzw.tum.de/coryneregnet/']
    doi = ['10.1038/s41597-020-0484-9']
    authors = [['Mariana Teixeira Dornelles Parise', 'Doglas Parise', 'Rodrigo Bentes Kato',
               'Josch Konstantin Pauling', 'Andreas Tauch', 'Vasco Ariston de Carvalho Azevedo', 'Jan Baumbach']]
    description = ['CoryneRegNet 7, the reference database and analysis platform for corynebacterial gene regulatory '
                   'networks']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(CoryneRegNetConnector,
                                source='coryneregnet',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self.create_connection(source='source', target='organism', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatorConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self.create_connection(source='source', target='regulator', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToGeneConnector(CoryneRegNetConnector,
                            source='coryneregnet',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self.create_connection(source='source', target='gene', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToTFBSConnector(CoryneRegNetConnector,
                            source='coryneregnet',
                            version='0.0.0',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):

    def connect(self):
        df = self.create_connection(source='source', target='tfbs', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(CoryneRegNetConnector,
                                             source='coryneregnet',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self.create_connection(source='source', target='rin', cardinality='one_to_many')
        self.stack_connections(df)
