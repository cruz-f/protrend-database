from typing import List

from protrend.model.model import Organism
from protrend.model.node import Node


class LoaderSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_files: List[str] = []
    default_node: Node = Node

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 node: Node = None,
                 files: List[str] = None):

        if not files:
            files = []

        self._source = source
        self._version = version
        self._node = node
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

    @property
    def node(self) -> Node:
        if not self._node:
            return self.default_node
