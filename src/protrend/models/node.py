from typing import Union, Type, TYPE_CHECKING

from neomodel import (UniqueIdProperty, RelationshipTo, DateTimeProperty)

if TYPE_CHECKING:
    from protrend.models.version import VersionNode


class Node:
    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    def __init_subclass__(cls,
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
    def cls_name(cls):
        return cls.__name__.lower()

    def connect_version(self, version: Union[None, VersionNode]):

        if version is None:
            return

        if hasattr(self, 'version'):
            self.version.connect(version)

        if hasattr(version, self.cls_name()):
            relationship = getattr(version, self.cls_name())
            relationship.connect(self)
