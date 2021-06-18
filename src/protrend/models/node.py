from typing import Union, Type, TYPE_CHECKING, Any

from neomodel import (UniqueIdProperty, RelationshipTo, DateTimeProperty, StructuredNode)

if TYPE_CHECKING:
    from protrend.models.version import VersionNode


class Node:
    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    def __init_subclass__(cls: Type[StructuredNode],
                          version: Type['VersionNode'] = None,
                          to_: bool = False,
                          from_: bool = False,
                          **kwargs):

        if version:

            if to_:
                setattr(cls, 'version', RelationshipTo(version, 'VERSIONING'))

            if from_:
                version.register(cls)

    @classmethod
    def cls_name(cls: Type[StructuredNode]):
        return cls.__name__.lower()

    @classmethod
    def properties(cls: Type[StructuredNode]):

        return cls.__all_properties__

    @classmethod
    def relationships(cls: Type[StructuredNode]):

        return cls.__all_relationships__

    @classmethod
    def from_item(cls: Type[StructuredNode], item: dict, save: bool = True) -> StructuredNode:
        properties = {attr: value for attr, value in item.items()
                      if attr in cls.properties()}

        instance = cls(**properties)

        if save:
            return instance.save()

        return instance

    @classmethod
    def get_by_version(cls: Type[StructuredNode], attr: str, value: Any, version: VersionNode):

        node_set = cls.nodes.filter(**{attr: value})

        for node in node_set:

            if node.version.is_connected(version):
                return node

    def update(self: StructuredNode, item: dict, save: bool = True):

        properties = {attr: value for attr, value in item.items()
                      if attr in self.properties()}

        for attr, value in item.items():

            if attr in self.properties():

                setattr(self, attr, value)

        if save:
            return self.save()

        return self

    def connect_version(self: StructuredNode, version: VersionNode):

        if hasattr(self, 'version'):
            self.version.connect(version)

        if hasattr(version, self.cls_name()):
            relationship = getattr(version, self.cls_name())
            relationship.connect(self)
