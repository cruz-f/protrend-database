from typing import TYPE_CHECKING, Any

from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode)

if TYPE_CHECKING:
    from protrend.models.version import Version


class Node(StructuredNode):

    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    __abstract_node__ = True

    @classmethod
    def cls_name(cls, lower: bool = False):

        if lower:
            return cls.__name__.lower()

        return cls.__name__

    @classmethod
    def properties(cls):

        return dict(cls.__all_properties__)

    @classmethod
    def relationships(cls):

        return dict(cls.__all_relationships__)

    @classmethod
    def from_item(cls, item: dict, save: bool = True) -> StructuredNode:
        properties = {attr: value for attr, value in item.items()
                      if attr in cls.properties()}

        instance = cls(**properties)

        if save:
            return instance.save()

        return instance

    @classmethod
    def get_by_version(cls, attr: str, value: Any, version: 'Version'):

        node_set = cls.nodes.filter(**{attr: value})

        for node in node_set:

            if node.version.is_connected(version):
                return node

    def update(self, item: dict, save: bool = True):

        properties = {attr: value for attr, value in item.items()
                      if attr in self.properties()}

        for attr, value in item.items():

            if attr in self.properties():

                setattr(self, attr, value)

        if save:
            return self.save()

        return self

    def connect_to_version(self, version: 'Version'):

        if hasattr(self, 'version'):
            self.version.connect(version)

        rel_name = self.cls_name(lower=True)

        if hasattr(version, rel_name):
            relationship = getattr(version, rel_name)
            relationship.connect(self)
