import os
from configparser import ConfigParser
from pathlib import Path
from typing import Union

from .singleton import singleton


@singleton
class Settings:

    def __init__(self):
        self._source = Path(os.path.dirname(__file__)).parent

        config = ConfigParser()
        config.read(self._source.joinpath('etl.conf'))

        self._working_directory = Path(config.get('etl-configuration', 'working_directory'))

        self._request_sleep = float(config.get('etl-configuration', 'request_sleep'))
        self._request_timeout = float(config.get('etl-configuration', 'request_timeout'))
        self._request_retries = int(config.get('etl-configuration', 'request_retries'))

        self._db_user_name = str(config.get('db-configuration', 'user_name'))
        self._db_password = str(config.get('db-configuration', 'password'))
        self._db_ip = str(config.get('db-configuration', 'ip'))
        self._db_port = str(config.get('db-configuration', 'port'))
        self._lasagna_url = str(config.get('db-configuration', 'lasagna_url'))

        self._started = False

    @property
    def source(self):
        return self._source

    @property
    def working_directory(self):
        return self._working_directory

    @property
    def request_sleep(self):
        return self._request_sleep

    @property
    def request_timeout(self):
        return self._request_timeout

    @property
    def request_retries(self):
        return self._request_retries

    @property
    def db_user_name(self):
        return self._db_user_name

    @property
    def db_password(self):
        return self._db_password

    @property
    def db_ip(self):
        return self._db_ip

    @property
    def db_port(self):
        return self._db_port

    @property
    def lasagna_url(self):
        return self._lasagna_url

    @property
    def extract(self):
        return self.source.joinpath('extract')

    @property
    def transform(self):
        return self.source.joinpath('transform')

    @property
    def load(self):
        return self.source.joinpath('load')

    @property
    def log(self):
        return self.source.joinpath('log')

    @property
    def data_lake(self):
        return self.working_directory.joinpath('data_lake')

    @property
    def bioapi_cache(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache')

    @property
    def entrez(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'entrez')

    @property
    def entrez_search(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'entrez', 'search')

    @property
    def entrez_summary(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'entrez', 'summary')

    @property
    def entrez_fetch(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'entrez', 'fetch')

    @property
    def kegg(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'kegg')

    @property
    def uniprot(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'uniprot')

    @property
    def uniprot_record(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'uniprot', 'record')

    @property
    def uniprot_query(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'uniprot', 'query')

    @property
    def uniprot_mapping(self):
        return self.working_directory.joinpath('data_lake', 'bioapi_cache', 'uniprot', 'mapping')

    @property
    def log_conf_file(self):
        return self.source.joinpath('log', 'log.conf')

    @property
    def log_working_directory(self):
        return self.working_directory.joinpath('log')

    def start_settings(self,
                       working_directory: Union[str, Path] = None,
                       request_sleep: float = None,
                       request_timeout: float = None,
                       request_retries: int = None):

        if working_directory:

            if isinstance(working_directory, str):
                working_directory = Path(working_directory)

            self._working_directory = working_directory

        if request_sleep is not None:
            self._request_sleep = request_sleep

        if request_timeout is not None:
            self._request_timeout = request_timeout

        if request_retries is not None:
            self._request_retries = request_retries

        self._started = True


Settings: Settings
