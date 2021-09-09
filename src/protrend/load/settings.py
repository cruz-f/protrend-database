from typing import List


class LoaderSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_files: List[str] = []

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 files: List[str] = None):

        if not files:
            files = []

        self._source = source
        self._version = version
        self._files = files

    @property
    def source(self) -> str:
        if not self._source:
            return self.default_source

    @property
    def version(self) -> str:
        if not self._version:
            return self.default_version

    @property
    def files(self) -> List[str]:
        if not self._files:
            return self.default_files
