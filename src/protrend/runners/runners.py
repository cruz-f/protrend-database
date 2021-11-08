from pathlib import Path

from protrend.utils import Settings, log_file_from_name
from protrend.log import ProtrendLogger
from protrend.pipeline import Pipeline
from protrend.utils import NeoDatabase


def run_database(install_labels: bool = False,
                 clear_constraints: bool = False,
                 clear_indexes: bool = False,
                 user_name: str = 'neo4j',
                 password: str = 'protrend',
                 ip: str = 'localhost',
                 port: str = '7687'):
    ProtrendLogger.log.info(f'Starting database {user_name} with {ip}:{port}')
    neo_db = NeoDatabase(user_name=user_name, password=password, ip=ip, port=port)
    neo_db.connect()

    if install_labels:
        ProtrendLogger.log.info(f'Installing new database labels')
        neo_db.auto_install_labels()

    if clear_constraints or clear_indexes:
        ProtrendLogger.log.info(f'Clearing old database constraints or indexes')
        neo_db.clear_db(clear_constraints=clear_constraints, clear_indexes=clear_indexes)

    return


def run_logger(name: str):
    log_file = log_file_from_name(name)
    ProtrendLogger.log_file = log_file
    ProtrendLogger.start_logger()
    return


def run_pipeline(source: str,
                 version: str,
                 extract: bool = False,
                 transform: bool = True,
                 connect: bool = True,
                 load: bool = True):
    if extract and transform and connect and load:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for ETCL')

        pipeline = Pipeline.for_etcl(source=source, version=version)
        pipeline.extract()
        pipeline.transform()
        pipeline.connect()
        pipeline.load()
        return

    if transform and connect and load:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for TCL')

        pipeline = Pipeline.for_tcl(source=source, version=version)
        pipeline.transform()
        pipeline.connect()
        pipeline.load()
        return

    if transform and connect:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for TC')

        pipeline = Pipeline.for_tc(source=source, version=version)
        pipeline.transform()
        pipeline.connect()
        return

    if extract:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for extraction')

        pipeline = Pipeline.for_extraction(source=source, version=version)
        pipeline.extract()
        return

    if transform:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for transformation')

        pipeline = Pipeline.for_transformation(source=source, version=version)
        pipeline.transform()
        return

    if connect:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for connection')

        pipeline = Pipeline.for_connection(source=source, version=version)
        pipeline.connect()
        return

    if load:
        ProtrendLogger.log.info(f'Starting pipeline with {source} and {version} for loading')

        pipeline = Pipeline.for_loading(source=source, version=version)
        pipeline.load()
        return


if __name__ == "__main__":
    # -------------------------------------------------------------
    # SET WORKING DIRECTORY FOR STAGING AREA AND DATA LAKE PATH
    # -------------------------------------------------------------
    working_directory = Path(r'C:\Users\BiSBII\Desktop\protrend')
    Settings.start_settings(working_directory=working_directory)

    run_logger('tcl_logger')
    run_database(install_labels=True, clear_constraints=True, clear_indexes=True)

    # ----------------------------------------------------
    # Abasy
    # ----------------------------------------------------
    # run_logger('abasy_logger')
    run_pipeline(source='abasy', version='0.0.0')

    # ----------------------------------------------------
    # CollecTF
    # ----------------------------------------------------
    # run_logger('collectf_logger')
    run_pipeline(source='collectf', version='0.0.1')

    # ----------------------------------------------------
    # CoryneRegNet
    # ----------------------------------------------------
    # run_logger('coryneregnet_logger')
    run_pipeline(source='coryneregnet', version='0.0.0')

    # ----------------------------------------------------
    # DBTBS
    # ----------------------------------------------------
    # run_logger('dbtbs_logger')
    run_pipeline(source='dbtbs', version='0.0.3')

    # ----------------------------------------------------
    # Literature
    # ----------------------------------------------------
    # run_logger('literature_logger')
    run_pipeline(source='literature', version='0.0.0')

    # ----------------------------------------------------
    # RegPrecise
    # ----------------------------------------------------
    # run_logger('regprecise_logger')
    run_pipeline(source='regprecise', version='0.0.0')

    # ----------------------------------------------------
    # RegulonDB
    # ----------------------------------------------------
    # run_logger('regulondb_logger')
    run_pipeline(source='regulondb', version='0.0.0')
