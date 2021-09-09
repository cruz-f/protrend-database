import sys
from pathlib import Path
from typing import Tuple, Dict

import pandas as pd

from protrend.io.json import read_json_frame
from protrend.load import RegPreciseLoader
from protrend.model.node import Node
from protrend.runners import Director
from protrend.transform.regprecise import *
from protrend.transform.regprecise.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils import NeoDatabase, ROOT_PATH, DATA_LAKE_PATH, STAGING_AREA_PATH

src_path = ROOT_PATH.parent
sys.path.insert(0, str(src_path))
from protrend.runners import run_spider


def extract_runner(spider: str = 'regprecise', staging_area: Path = STAGING_AREA_PATH, version: str = '0.0.0'):
    return run_spider(spider=spider, staging_area=staging_area, version=version)


def transform_runner(transform: bool = True,
                     connect: bool = True,
                     install_labels: bool = False,
                     clear_constraints: bool = False,
                     clear_indexes: bool = False) -> Tuple[Director, Dict[str, pd.DataFrame]]:
    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    if install_labels:
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    transformers = [
        EffectorTransformer(),
        GeneTransformer(),
        OperonTransformer(),
        OrganismTransformer(),
        PathwayTransformer(),
        PublicationTransformer(),
        RegulatorTransformer(),
        RegulatoryFamilyTransformer(),
        RegulatoryInteractionTransformer(),
        SourceTransformer(),
        TFBSTransformer(),
    ]
    connectors = [
        EffectorToOrganismConnector(), EffectorToRegulatorConnector(), EffectorToSourceConnector(),
        GeneToOrganismConnector(), GeneToSourceConnector(), GeneToRegulatorConnector(), GeneToTFBSConnector(),
        OperonToGeneConnector(), OperonToOrganismConnector(), OperonToRegulatorConnector(), OperonToSourceConnector(),
        OperonToTFBSConnector(),
        OrganismToSourceConnector(),
        PathwayToGeneConnector(), PathwayToRegulatorConnector(), PathwayToSourceConnector(),
        RegulatorToOrganismConnector(), RegulatorToSourceConnector(),
        RegulatoryFamilyToPublicationConnector(), RegulatoryFamilyToRegulatorConnector(),
        RegulatoryFamilyToSourceConnector(),
        RegulatoryInteractionToEffectorConnector(), RegulatoryInteractionToGeneConnector(),
        RegulatoryInteractionToOperonConnector(), RegulatoryInteractionToOrganismConnector(),
        RegulatoryInteractionToRegulatorConnector(), RegulatoryInteractionToSourceConnector(),
        RegulatoryInteractionToTFBSConnector(),
        TFBSToRegulatorConnector(), TFBSToOrganismConnector(), TFBSToSourceConnector(),
    ]

    director = Director(transformers=transformers,
                        connectors=connectors)

    if transform:
        director.transform()

    if connect:
        director.connect()

    regprecise_data_lake = DATA_LAKE_PATH.joinpath('regprecise', '0.0.0')
    data_lake_files = regprecise_data_lake.glob('*.json')
    data_lake = {fp.stem: read_json_frame(fp)
                 for fp in data_lake_files}

    return director, data_lake


def load_runner(install_labels: bool = False,
                clear_constraints: bool = False,
                clear_indexes: bool = False) -> Tuple[Director, Dict[str, pd.DataFrame]]:
    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    if install_labels:
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    loaders = [RegPreciseLoader()]
    director = Director(loaders=loaders)
    director.load()

    database = {node_name: node.node_to_df()
                for node_name, node in Node.node_register.items()}

    return director, database


if __name__ == "__main__":
    # ----------------------------------------------------
    # EXTRACT
    # ----------------------------------------------------
    # extract_runner()

    # ----------------------------------------------------
    # TRANSFORM
    # ----------------------------------------------------
    reg_transformer, reg_data_lake = transform_runner()

    # ----------------------------------------------------
    # LOAD
    # ----------------------------------------------------
    reg_loader, reg_database = load_runner()
