from protrend.model import Source, Organism, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform import BaseSourceTransformer
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.utils import SetList


class SourceTransformer(CoryneRegNetTransformer, BaseSourceTransformer,
                        source='coryneregnet',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['coryneregnet']
    type = ['database']
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
    default_connect_stack = {'source': 'integrated_source.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self.create_connection(source='source', target='organism', cardinality='one_to_many')
        self.stack_json(df)


class SourceToRegulatorConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self.create_connection(source='source', target='regulator', cardinality='one_to_many')
        self.stack_json(df)


class SourceToGeneConnector(CoryneRegNetConnector,
                            source='coryneregnet',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='source', target='gene', cardinality='one_to_many')
        self.stack_json(df)


class SourceToTFBSConnector(CoryneRegNetConnector,
                            source='coryneregnet',
                            version='0.0.0',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='source', target='tfbs', cardinality='one_to_many')
        self.stack_json(df)


class SourceToRegulatoryInteractionConnector(CoryneRegNetConnector,
                                             source='coryneregnet',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='source', target='rin', cardinality='one_to_many')
        self.stack_json(df)
