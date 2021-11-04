from protrend.load import AbasyLoader
from protrend.runners import Director
from protrend.utils import NeoDatabase


def transform_load_abasy(install_labels, clear_constraints, clear_indexes):

    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    if install_labels:
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    transformers = [
        GeneTransformer(),
        OperonTransformer(),
        OrganismTransformer(),
        RegulatorTransformer(),
        RegulatoryInteractionTransformer(),
        SourceTransformer(),
    ]
    connectors = [
        OperonToGeneConnector(),

        GeneToOrganismConnector(),
        OperonToOrganismConnector(),
        RegulatorToOrganismConnector(),
        RegulatoryInteractionToOrganismConnector(),

        RegulatorToGeneConnector(),
        RegulatorToOperonConnector(),

        RegulatoryInteractionToGeneConnector(),
        RegulatoryInteractionToOperonConnector(),
        RegulatoryInteractionToRegulatorConnector(),

        GeneToSourceConnector(),
        OperonToSourceConnector(),
        OrganismToSourceConnector(),
        RegulatorToSourceConnector(),
        RegulatoryInteractionToSourceConnector(),
    ]

    loaders = [AbasyLoader()]

    director = Director(transformers=transformers,
                        connectors=connectors,
                        loaders=loaders)

    director.transform()
    director.connect()
    director.load()