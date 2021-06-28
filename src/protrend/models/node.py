from typing import TYPE_CHECKING, Any, Union, Callable

import pandas as pd
from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode)

if TYPE_CHECKING:
    from protrend.models.version import Version


class Node(StructuredNode):

    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    __abstract_node__ = True

    @classmethod
    def cls_name(cls, transform: Callable):

        if transform:
            return transform(cls.__name__)

        return cls.__name__

    @classmethod
    def cls_properties(cls):

        return dict(cls.__all_properties__)

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
    def cls_to_dict(cls):

        res = {key: [] for key in cls.cls_keys()}

        for node in cls.nodes.all():
            for key, val in node.properties.items():
                if key in res:
                    res[key].append(val)

        return res

    @classmethod
    def cls_to_df(cls):
        return pd.DataFrame.from_dict(cls.cls_to_dict())

    @classmethod
    def from_item(cls, item: dict, save: bool = True) -> StructuredNode:
        properties = {attr: value for attr, value in item.items()
                      if attr in cls.cls_properties()}

        instance = cls(**properties)

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
    def properties(self):
        return {key: val for key, val in self.__properties__.items()
                if val is not None}

    def update(self, item: dict, save: bool = True) -> 'Node':

        properties = {attr: value for attr, value in item.items()
                      if attr in self.cls_properties()}

        for attr, value in item.items():

            if attr in self.cls_properties():

                setattr(self, attr, value)

        if save:
            return self.save()

        return self

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
