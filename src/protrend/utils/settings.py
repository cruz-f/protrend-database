import os
from pathlib import Path


class Settings:
    REQUEST_SLEEP: float = 0.25
    REQUEST_TIMEOUT: int = 30
    REQUEST_RETRIES: int = 3

    ROOT_PATH: Path = Path(os.path.dirname(__file__)).parent
    EXTRACT_PATH: Path = ROOT_PATH.joinpath('extract')
    TRANSFORM_PATH: Path = ROOT_PATH.joinpath('transform')
    LOAD_PATH: Path = ROOT_PATH.joinpath('load')
    STAGING_AREA_PATH: Path = ROOT_PATH.joinpath('staging_area')
    DATA_LAKE_PATH: Path = ROOT_PATH.joinpath('data_lake')
    DATA_LAKE_BIOAPI_PATH: Path = ROOT_PATH.joinpath('data_lake', 'bioapi_cache')
