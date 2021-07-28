from collections import defaultdict
from datetime import datetime
from typing import List, Dict, Any, Union

import pandas as pd
import pytz
from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode, StringProperty)


class Node(StructuredNode):
    __abstract_node__ = True

    uid = UniqueIdProperty()
    protrend_id = StringProperty(required=True)
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    identifying_property = 'protrend_id'
    header = 'PRT'
    entity = 'PRT'

    # -------------------------------------
    # Class attributes
    # -------------------------------------
    @classmethod
    def node_name(cls):
        return cls.__name__.lower()

    @classmethod
    def node_properties(cls) -> dict:
        return dict(cls.__all_properties__)

    @classmethod
    def node_relationships(cls) -> dict:
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
    def latest_identifier(cls):
        protrend_identifiers = [node.identifier for node in cls.nodes.all()]
        if protrend_identifiers:
            sorted_identifiers = sorted(protrend_identifiers, key=protrend_id_decoder, reverse=True)
            return sorted_identifiers[0]

        return protrend_id_encoder(cls.header, cls.entity, 0)

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
    def node_to_dict(cls, to: str = 'dict') -> dict:

        if to == 'dict':
            res = defaultdict(list)

            for node in cls.nodes.all():

                node_properties = node.properties

                for key in cls.node_keys():
                    val = node_properties.get(key, None)
                    res[key].append(val)

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
                           if key in node_keys and val is not None}

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

                    if key in node_keys and val is not None:
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


def protrend_id_encoder(header: str, entity: str, integer: Union[str, int]):

    integer = int(integer)

    return f'{header}.{entity}.{integer:7d}'


def protrend_id_decoder(protrend_id: str):

    prt, entity, integer = protrend_id.split('.')

    return int(integer)