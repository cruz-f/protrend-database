import os
from pathlib import Path

REQUEST_SLEEP: float = 0.25
REQUEST_TIMEOUT: int = 30
REQUEST_RETRIES: int = 3

ROOT_PATH: Path = Path(os.path.dirname(__file__)).parent
EXTRACT_PATH: Path = ROOT_PATH.joinpath('extract')
STAGING_AREA_PATH: Path = ROOT_PATH.joinpath('extract', 'staging_area')
TRANSFORM_PATH: Path = ROOT_PATH.joinpath('transform')
DATA_LAKE_PATH: Path = ROOT_PATH.joinpath('transform', 'data_lake')
DATA_LAKE_BIOAPI_PATH: Path = ROOT_PATH.joinpath('transform', 'data_lake', 'bioapi_cache')