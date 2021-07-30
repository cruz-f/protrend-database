from typing import Dict, Tuple

from protrend.model.node import Node


class TransformerSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_files: Dict[str, str] = {}
    default_node: Node = Node
    default_node_factors: Tuple[str] = ()

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 node: Node = None,
                 node_factors: Tuple[str] = None,
                 **files: Dict[str, str]):

        self._source = source
        self._version = version
        self._node = node
        self._node_factors = node_factors
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
    def files(self) -> Dict[str, str]:
        if not self._files:
            return self.default_files

    @property
    def node(self) -> Node:
        if not self._node:
            return self.default_node

    @property
    def node_factors(self) -> Tuple[str]:
        if not self._node_factors:
            return self.default_node_factors
