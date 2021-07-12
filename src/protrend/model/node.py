from datetime import datetime
from typing import Callable, List

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

    # -------------------------------------
    # Class attributes
    # -------------------------------------
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

    # -------------------------------------
    # Class dynamic attributes/methods
    # -------------------------------------
    @classmethod
    def cls_to_dict(cls, properties: List[str] = None) -> dict:

        if not properties:
            res = {key: [] for key in cls.cls_keys()}

        else:
            cls_properties = list(cls.cls_keys())
            res = {key: [] for key in properties
                   if key in cls_properties}

        for node in cls.nodes.all():
            for key, val in node.properties.items():
                if key in res:
                    res[key].append(val)

        return res

    @classmethod
    def cls_to_df(cls, properties: List[str] = None) -> pd.DataFrame:
        return pd.DataFrame.from_dict(cls.cls_to_dict(properties))

    # -------------------------------------
    # Class extensions
    # -------------------------------------
    @classmethod
    def create_or_update(cls, *props: dict, **kwargs):

        for prop in props:
            prop['updated'] = datetime.utcnow().replace(tzinfo=pytz.utc)

        return super(Node, cls).create_or_update(*props, **kwargs)

    # -------------------------------------
    # Instance dynamic attributes/methods
    # -------------------------------------
    @property
    def identifier(self):
        return self.properties[self.identifying_property]

    @property
    def properties(self) -> dict:
        return {key: val for key, val in self.__properties__.items()
                if val is not None}

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
