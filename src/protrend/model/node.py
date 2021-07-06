from typing import TYPE_CHECKING, Any, Union, Callable

import pandas as pd
from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode)

if TYPE_CHECKING:
    from protrend.model.version import Version


class Node(StructuredNode):

    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    property_as_id = 'uid'

    __abstract_node__ = True

    def __init__(self, *args, **kwargs):
        self._relationships = {}
        super(Node, self).__init__(*args, **kwargs)

    @classmethod
    def cls_name(cls, transform: Callable = None):

        if transform:
            return transform(cls.__name__)

        return cls.__name__

    @classmethod
    def cls_properties(cls) -> dict:

        return dict(cls.__all_properties__)

    @classmethod
    def cls_relationships(cls) -> dict:

        return dict(cls.__all_relationships__)

    @classmethod
    def cls_keys(cls):
        return cls.cls_properties().keys()

    @classmethod
    def cls_values(cls):
        return cls.cls_properties().values()

    @classmethod
    def cls_items(cls):
        return cls.cls_properties().items()

    @classmethod
    def cls_to_dict(cls) -> dict:

        res = {key: [] for key in cls.cls_keys()}

        for node in cls.nodes.all():
            for key, val in node.properties.items():
                if key in res:
                    res[key].append(val)

        return res

    @classmethod
    def cls_to_df(cls) -> pd.DataFrame:
        return pd.DataFrame.from_dict(cls.cls_to_dict())

    @classmethod
    def from_item(cls, item: dict, save: bool = True) -> 'Node':

        cls_properties = cls.cls_properties()
        cls_relationships = cls.cls_relationships()

        properties = {}
        relationships = {}

        for key, val in item.items():
            if key in cls_properties:
                properties[key] = val

            elif key in cls_relationships:
                relationships[key] = val

        instance = cls(**properties)

        instance.relationships = relationships

        if save:
            return instance.save()

        return instance

    @classmethod
    def get_by_version(cls, attr: str, value: Any, version: 'Version') -> Union['Node', None]:

        node_set = cls.nodes.filter(**{attr: value})

        for node in node_set:

            if node.version.is_connected(version):
                return node

    @property
    def identifier(self):
        return self.properties[self.property_as_id]

    @property
    def properties(self) -> dict:
        return {key: val for key, val in self.__properties__.items()
                if val is not None}

    @property
    def relationships(self) -> dict:
        return self._relationships.copy()

    @relationships.setter
    def relationships(self, value: Union[dict, tuple]):

        if isinstance(value, dict):

            self._relationships = {}

            for key, val in value.items():
                if isinstance(val, (tuple, list, set)):
                    self._relationships[key] = val

                else:
                    self._relationships[key] = [val]

        elif isinstance(value, tuple):
            new_rel = {value[0]: value[1]}
            self._relationships.update(new_rel)

    def connect_to_version(self, version: 'Version'):

        if hasattr(self, 'version'):
            self.version.connect(version)

        rel_name = self.cls_name(str.lower)

        if hasattr(version, rel_name):
            relationship = getattr(version, rel_name)
            relationship.connect(self)

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

    def relationship_to_df(self, relationship: str):
        data = [[self.identifier, val] for val in self.relationships.get(relationship, [])]
        return pd.DataFrame(data, columns=['start', 'end'])
