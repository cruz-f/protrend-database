import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model import Source, Organism, Regulator, Gene, RegulatoryInteraction
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.abasy.organism import OrganismTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.transform.abasy.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils import SetList


class SourceTransformer(AbasyTransformer,
                        source='abasy',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = 'abasy'
    type = 'database'
    url = 'https://abasy.ccg.unam.mx'
    doi = '10.1016/j.csbj.2020.05.015'
    authors = ['Juan M. Escorcia-Rodríguez, AndreasTauch, Julio A. Freyre-Gonzáleza']
    description = 'Abasy Atlas v2.2: The most comprehensive and up-to-date inventory of meta-curated, historical, bacterial regulatory networks, their completeness and system-level characterization'

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self):
        db = dict(name=[self.name],
                  type=[self.type],
                  url=[self.url],
                  doi=[self.doi],
                  authors=[self.authors],
                  description=[self.description])

        df = pd.DataFrame(db, index=[0])

        self.stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(AbasyConnector,
                                source='abasy',
                                version='0.0.0',
                                from_node=Organism,
                                to_node=Source,
                                register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'source': 'integrated_source.json'}

    def connect(self):
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = organism['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToSourceConnector(AbasyConnector,
                                 source='abasy',
                                 version='0.0.0',
                                 from_node=Regulator,
                                 to_node=Source,
                                 register=True):
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'source': 'integrated_source.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToSourceConnector(AbasyConnector,
                            source='abasy',
                            version='0.0.0',
                            from_node=Gene,
                            to_node=Source,
                            register=True):
    default_connect_stack = {'gene': 'integrated_gene.json', 'source': 'integrated_source.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToSourceConnector(AbasyConnector,
                                             source='abasy',
                                             version='0.0.0',
                                             from_node=RegulatoryInteraction,
                                             to_node=Source,
                                             register=True):

    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
