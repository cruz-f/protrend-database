import os
from pathlib import Path
from typing import Union

from .singleton import singleton


@singleton
class Settings:

    def __init__(self):
        self._source = Path(os.path.dirname(__file__)).parent
        self._working_directory = Path(os.path.dirname(__file__)).parent

        self._request_sleep = 0.25
        self._request_timeout = 30
        self._request_retries = 3

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
    def staging_area(self):
        return self.working_directory.joinpath('staging_area')

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

