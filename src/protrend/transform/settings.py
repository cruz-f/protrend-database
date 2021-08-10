from typing import Dict, Tuple, Type

from protrend.model.node import Node


class TransformerSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_node: Type[Node] = Node
    default_node_factors: Tuple[str] = ()
    default_transform: Dict[str, str] = {}
    default_order: int = 0

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 node: Type[Node] = None,
                 node_factors: Tuple[str] = None,
                 transform: Dict[str, str] = None,
                 order: int = None):

        if not transform:
            transform = {}

        self._source = source
        self._version = version
        self._node = node
        self._node_factors = node_factors
        self._transform = transform
        self._order = order

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.node.node_name()} - {self.node_factors}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'TransformerSettings'):

        other_values = {other.__class__.__name__, other.source, other.version, other.node.node_name(),
                        other.node_factors}

        if self.__class__.__name__ not in other_values:
            return False

        if self.source not in other_values:
            return False

        if self.version not in other_values:
            return False

        if self.node.node_name() not in other_values:
            return False

        if self.node_factors not in other_values:
            return False

        return True

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
    def node(self) -> Type[Node]:
        if not self._node:
            return self.default_node

        return self._node

    @property
    def node_factors(self) -> Tuple[str]:
        if not self._node_factors:
            return self.default_node_factors

        return self._node_factors

    @property
    def transform(self) -> Dict[str, str]:

        if not self._transform:
            return self.default_transform

        return self._transform

    @property
    def order(self) -> int:
        if self._order is None:
            return self.default_order

        return self._order


class ConnectorSettings:
    default_source: str = ''
    default_version: str = '0.0.0'
    default_from_node: Type[Node] = Node
    default_to_node: Type[Node] = Node
    default_connect: Dict[str, str] = {}

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 from_node: Type[Node] = None,
                 to_node: Type[Node] = None,
                 connect: Dict[str, str] = None):

        if not connect:
            connect = {}

        self._source = source
        self._version = version
        self._from_node = from_node
        self._to_node = to_node
        self._connect = connect

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.from_node.node_name()} - {self.to_node.node_name()}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'ConnectorSettings'):

        other_values = {other.__class__.__name__, other.source, other.version, other.from_node.node_name(),
                        other.to_node.node_name()}

        if self.__class__.__name__ not in other_values:
            return False

        if self.source not in other_values:
            return False

        if self.version not in other_values:
            return False

        if self.from_node.node_name() not in other_values:
            return False

        if self.to_node.node_name() not in other_values:
            return False

        return True

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
    def from_node(self) -> Type[Node]:
        if not self._from_node:
            return self.default_from_node

        return self._from_node

    @property
    def to_node(self) -> Type[Node]:
        if not self._to_node:
            return self.default_to_node

        return self._to_node

    @property
    def connect(self) -> Dict[str, str]:
        if not self._connect:
            return self.default_connect

        return self._connect
