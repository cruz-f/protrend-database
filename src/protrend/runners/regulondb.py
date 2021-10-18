from typing import Tuple, Dict

import pandas as pd

from protrend.io.json import read_json_frame
from protrend.load import RegulondbLoader
from protrend.log import ProtrendLogger
from protrend.model.node import Node
from protrend.runners import Director
from protrend.transform.regulondb import *
from protrend.utils import NeoDatabase, DATA_LAKE_PATH
from protrend.utils.miscellaneous import log_file_from_name


def transform_runner(transform: bool = True,
                     connect: bool = True,
                     install_labels: bool = False,
                     clear_constraints: bool = False,
                     clear_indexes: bool = False) -> Tuple[Director, Dict[str, pd.DataFrame]]:
    ProtrendLogger.log.info(f'Starting transform runner with transform: {transform}, connect: {connect}, '
                            f'install labels: {install_labels}, clear constraints: {clear_constraints}, '
                            f'clear indexes: {clear_indexes}')

    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    if install_labels:
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    transformers = [
        EffectorTransformer(),
        EvidenceTransformer(),
        GeneTransformer(),
        OperonTransformer(),
        OrganismTransformer(),
        PromoterTransformer(),
        PublicationTransformer(),
        RegulatorTransformer(),
        RegulatoryFamilyTransformer(),
        RegulatoryInteractionTransformer(),
        SourceTransformer(),
        TFBSTransformer(),
    ]
    connectors = [
        EvidenceToGeneConnector(),
        EvidenceToOperonConnector(),
        EvidenceToPromoterConnector(),
        EvidenceToRegulatorConnector(),
        EvidenceToRegulatoryInteractionConnector(),
        EvidenceToTFBSConnector(),

        GeneToPromoterConnector(),
        GeneToTFBSConnector(),
        OperonToGeneConnector(),
        OperonToPromoterConnector(),
        OperonToTFBSConnector(),

        EffectorToOrganismConnector(),
        GeneToOrganismConnector(),
        OperonToOrganismConnector(),
        PromoterToOrganismConnector(),
        RegulatorToOrganismConnector(),
        RegulatoryInteractionToOrganismConnector(),
        TFBSToOrganismConnector(),

        PublicationToGeneConnector(),
        PublicationToOperonConnector(),
        PublicationToPromoterConnector(),
        PublicationToRegulatorConnector(),
        PublicationToTFBSConnector(),

        RegulatorToRegulatoryFamilyConnector(),

        RegulatorToEffectorConnector(),
        RegulatorToGeneConnector(),
        RegulatorToOperonConnector(),
        RegulatorToTFBSConnector(),
        RegulatoryInteractionToEffectorConnector(),
        RegulatoryInteractionToGeneConnector(),
        RegulatoryInteractionToOperonConnector(),
        RegulatoryInteractionToRegulatorConnector(),
        RegulatoryInteractionToTFBSConnector(),

        EffectorToSourceConnector(),
        GeneToSourceConnector(),
        OperonToSourceConnector(),
        OrganismToSourceConnector(),
        PromoterToSourceConnector(),
        RegulatorToSourceConnector(),
        RegulatoryFamilyToSourceConnector(),
        RegulatoryInteractionToSourceConnector(),
        TFBSToSourceConnector(),
    ]

    director = Director(transformers=transformers,
                        connectors=connectors)

    if transform:
        director.transform()

    if connect:
        director.connect()

    regulondb_data_lake = DATA_LAKE_PATH.joinpath('regulondb', '0.0.0')
    data_lake_files = regulondb_data_lake.glob('*.json')
    data_lake = {fp.stem: read_json_frame(fp)
                 for fp in data_lake_files}

    data_lake_info = [f'{key}: {df.shape}' for key, df in data_lake.items()]
    data_lake_info = ' '.join(data_lake_info)
    ProtrendLogger.log.info(f'Transform stats: {data_lake_info}')
    ProtrendLogger.log.info(f'Finished transform runner')

    return director, data_lake


def load_runner(install_labels: bool = False,
                clear_constraints: bool = False,
                clear_indexes: bool = False) -> Tuple[Director, Dict[str, pd.DataFrame]]:
    ProtrendLogger.log.info(f'Starting loader runner with transform: '
                            f'install labels: {install_labels}, clear constraints: {clear_constraints}, '
                            f'clear indexes: {clear_indexes}')

    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    if install_labels:
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    loaders = [RegulondbLoader()]
    director = Director(loaders=loaders)
    director.load()

    database = {node_name: node.node_to_df()
                for node_name, node in Node.node_register.items()}

    database_info = [f'{key}: {df.shape}' for key, df in database.items()]
    database_info = ' '.join(database_info)
    ProtrendLogger.log.info(f'Database stats: {database_info}')
    ProtrendLogger.log.info(f'Finished loader runner')

    return director, database


if __name__ == "__main__":
    # ----------------------------------------------------
    # EXTRACT
    # ----------------------------------------------------
    # extract_runner()

    # ----------------------------------------------------
    # LOGGER
    # ----------------------------------------------------
    log_file = log_file_from_name('regulondb_transform')
    ProtrendLogger.log_file = log_file
    ProtrendLogger.start_logger()

    # ----------------------------------------------------
    # TRANSFORM
    # ----------------------------------------------------
    reg_transformer, reg_data_lake = transform_runner(install_labels=True, clear_constraints=True,
                                                      clear_indexes=True)

    # ----------------------------------------------------
    # LOAD
    # ----------------------------------------------------
    reg_loader, reg_database = load_runner(install_labels=True, clear_constraints=True, clear_indexes=True)