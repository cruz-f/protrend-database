from typing import Union, Type, TYPE_CHECKING

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
    def from_item(cls: Type[StructuredNode], item: dict):
        pass

    def connect_version(self: StructuredNode, version: VersionNode):

        if hasattr(self, 'version'):
            self.version.connect(version)

        if hasattr(version, self.cls_name()):
            relationship = getattr(version, self.cls_name())
            relationship.connect(self)
