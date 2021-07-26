import os
from pathlib import Path

REQUEST_TIMEOUT: int = 30
REQUEST_RETRIES: int = 3

ROOT_PATH = Path(os.path.dirname(__file__)).parent
EXTRACT_PATH = ROOT_PATH.joinpath('extract')
STAGING_AREA_PATH: str = ROOT_PATH.joinpath('extract', 'staging_area')
TRANSFORM_PATH = ROOT_PATH.joinpath('transform')
DATA_LAKE_PATH = ROOT_PATH.joinpath('transform', 'data_lake')