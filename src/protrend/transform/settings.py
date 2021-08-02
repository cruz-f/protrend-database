from datetime import datetime
from typing import Dict, Tuple, List

import pytz

from protrend.model.node import Node


class RelationshipSettings:

    def __init__(self,
                 from_node: Node = None,
                 to_node: Node = None,
                 from_property: str = None,
                 to_property: str = None,
                 **properties):

        if not properties:
            properties = {}

        self.from_node = from_node
        self.to_node = to_node
        self.from_property = from_property
        self.to_property = to_property
        self.properties = properties

    def is_empty(self) -> bool:
        if self.from_node is None:
            return True

        return False

    @property
    def key(self):
        if self.is_empty():
            return f'connected_{datetime.utcnow().replace(tzinfo=pytz.utc)}'

        return f'connected_{self.from_node.node_name()}_{self.to_node.node_name()}'

    @property
    def tuple_key(self):
        if self.is_empty():
            return (f'relationship_{datetime.utcnow().replace(tzinfo=pytz.utc)}',
                    f'relationship_{datetime.utcnow().replace(tzinfo=pytz.utc)}')

        return self.from_node.node_name(), self.to_node.node_name()


class TransformerSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_files: Dict[str, str] = {}
    default_node: Node = Node
    default_node_factors: Tuple[str] = ()
    default_relationships: List[RelationshipSettings] = [RelationshipSettings()]

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 node: Node = None,
                 node_factors: Tuple[str] = None,
                 relationships: List[RelationshipSettings] = None,
                 files: Dict[str, str] = None):

        if not relationships:
            relationships = []

        if not files:
            files = {}

        self._source = source
        self._version = version
        self._node = node
        self._node_factors = node_factors
        self._relationships = relationships
        self._files = files

    @property
    def source(self) -> str:
        if not self._source:
            return self.default_source

        return self._source

    @property
    def version(self) -> str:
        if not self._version:
            return self.default_version

        return self._version

    @property
    def node(self) -> Node:
        if not self._node:
            return self.default_node

        return self._node

    @property
    def node_factors(self) -> Tuple[str]:
        if not self._node_factors:
            return self.default_node_factors

        return self._node_factors

    @property
    def files(self) -> Dict[str, str]:
        if not self._files:
            return self.default_files

        return self._files

    @property
    def relationships(self) -> List[RelationshipSettings]:

        if not self._relationships:
            return self.default_relationships

        return self._relationships
