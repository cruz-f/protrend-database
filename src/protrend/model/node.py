from collections import defaultdict
from datetime import datetime
from typing import List, Dict, Any, Union, Type, Tuple

import pandas as pd
import pytz
from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode, StringProperty, RelationshipManager)
from neomodel.relationship import RelationshipMeta, StructuredRel

from protrend.log.logger import Logger
from protrend.utils.miscellaneous import convert_to_snake_case, is_null


class Node(StructuredNode):
    __abstract_node__ = True

    uid = UniqueIdProperty()
    protrend_id = StringProperty(required=True)
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    identifying_property = 'protrend_id'
    header = 'PRT'
    entity = 'PRT'

    node_register: Dict[str, Type['Node']] = {}

    def __init_subclass__(cls, **kwargs):
        cls.node_register[cls.node_name()] = cls

    # -------------------------------------
    # Class attributes
    # -------------------------------------
    @classmethod
    def node_name(cls):
        return convert_to_snake_case(cls.__name__.lower())

    @classmethod
    def node_properties(cls) -> dict:
        return dict(cls.__all_properties__)

    @classmethod
    def node_relationships(cls) -> Dict[str, RelationshipManager]:
        return dict(cls.__all_relationships__)

    @classmethod
    def node_keys(cls):
        return cls.node_properties().keys()

    @classmethod
    def node_values(cls):
        return cls.node_properties().values()

    @classmethod
    def node_items(cls):
        return cls.node_properties().items()

    @classmethod
    def last_node(cls) -> Union['Node', None]:
        nodes = cls.nodes.all()
        if nodes:
            sorted_nodes = sorted(nodes, key=_sort_nodes, reverse=True)
            return sorted_nodes[0]

        return None

    # -------------------------------------
    # Class methods
    # -------------------------------------
    @classmethod
    def node_from_dict(cls,
                       *nodes: Dict[str, Any],
                       save: bool = False) -> List['Node']:

        structured_nodes = []
        node_keys = list(cls.node_keys())

        for node in nodes:

            node_kwargs = {key: val for key, val in node.items()
                           if key in node_keys and val is not None}

            if node_kwargs:

                structured_node = cls(**node_kwargs)

                if save:
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_update_from_dict(cls,
                              *nodes: Dict[str, Any],
                              save: bool = False) -> List['Node']:

        structured_nodes = []
        node_keys = [key for key in cls.node_keys() if key != cls.identifying_property]
        all_nodes = cls.node_to_dict(to='node')

        for node in nodes:

            identifier = node.get(cls.identifying_property, '')

            if identifier in all_nodes:

                structured_node = all_nodes[identifier]

                for key, val in node.items():

                    if key in node_keys and val is not None:
                        setattr(structured_node, key, val)

                if save:
                    structured_node.updated = datetime.utcnow().replace(tzinfo=pytz.utc)
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_to_dict(cls, to: str = 'dict') -> Union[Dict[str, Any], Dict[str, 'Node']]:

        if to == 'dict':
            res = defaultdict(list)

            for node in cls.nodes.all():

                node_properties = node.properties

                for key in cls.node_keys():
                    val = node_properties.get(key, None)
                    res[key].append(val)

            if not res:
                return {key: [] for key in cls.node_keys()}

            return res

        elif to == 'node':

            return {node.identifier: node for node in cls.nodes.all()}

        raise ValueError(f'Invalid output {to}')

    @classmethod
    def node_from_df(cls,
                     nodes: pd.DataFrame,
                     save: bool = False) -> List['Node']:

        structured_nodes = []
        node_keys = list(cls.node_keys())

        for _, node in nodes.iterrows():

            node_kwargs = {key: val for key, val in node.items()
                           if key in node_keys and not is_null(val)}

            if node_kwargs:

                structured_node = cls(**node_kwargs)

                if save:
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_update_from_df(cls,
                            nodes: pd.DataFrame,
                            save: bool = False) -> List['Node']:

        structured_nodes = []
        node_keys = [key for key in cls.node_keys() if key != cls.identifying_property]
        all_nodes = cls.node_to_dict(to='node')

        for _, node in nodes.iterrows():

            identifier = node.get(cls.identifying_property, '')

            if identifier in all_nodes:

                structured_node = all_nodes[identifier]

                for key, val in node.items():

                    if key in node_keys and not is_null(val):
                        setattr(structured_node, key, val)

                if save:
                    structured_node.updated = datetime.utcnow().replace(tzinfo=pytz.utc)
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_to_df(cls) -> pd.DataFrame:
        return pd.DataFrame.from_dict(cls.node_to_dict(to='dict'))

    # -------------------------------------
    # Instance methods
    # -------------------------------------
    @property
    def identifier(self):
        return self.properties[self.identifying_property]

    @property
    def properties(self) -> Dict[str, Any]:
        return {key: val for key, val in self.__properties__.items()
                if val is not None and key != 'id'}

    def keys(self):
        return self.properties.keys()

    def values(self):
        return self.properties.values()

    def items(self):
        return self.properties.items()

    def to_dict(self):
        return self.properties

    def to_series(self):
        return pd.Series(data=self.values(), index=self.keys())

    def to_df(self):
        data = [list(self.values())]
        return pd.DataFrame(data=data, columns=self.keys())


def _sort_nodes(node: Node):
    return protrend_id_decoder(node.protrend_id)


def protrend_id_encoder(header: str, entity: str, integer: Union[str, int]) -> str:
    integer = int(integer)

    return f'{header}.{entity}.{integer:07}'


def protrend_id_decoder(protrend_id: str) -> int:
    prt, entity, integer = protrend_id.split('.')

    return int(integer)


def get_node_by_name(name: str, default=None) -> Union[Type[Node], None]:
    return Node.node_register.get(name, default)


def _find_to_node(relationship: RelationshipManager) -> Type[Node]:
    return relationship.definition['node_class']


def get_nodes_relationships(from_node: Type[Node], to_node: Type[Node], default=None) -> Tuple[List[str], List[str]]:
    from_node_rels = from_node.node_relationships()
    from_node_matches = []
    for attr, relationship in from_node_rels.items():

        this_to_node = _find_to_node(relationship)
        if this_to_node.node_name() == to_node.node_name():
            from_node_matches.append(attr)

    to_node_rels = to_node.node_relationships()
    to_node_matches = []
    for attr, relationship in to_node_rels.items():

        this_from_node = _find_to_node(relationship)
        if this_from_node.node_name() == from_node.node_name():
            to_node_matches.append(attr)

    return from_node_matches, to_node_matches


def connect_nodes(from_node: Node, to_node: Node, relationship: str, kwargs: dict) -> bool:
    relationship: RelationshipManager = getattr(from_node, relationship, None)
    relationship_model: StructuredRel = relationship.definition['model']

    if kwargs:
        kwargs = {key: val for key, val in kwargs.items()
                  if hasattr(relationship_model, key)}

    else:
        kwargs = {}

    return relationship.connect(to_node, kwargs)
