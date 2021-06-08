from typing import Type, TYPE_CHECKING, Union, Dict

from neomodel import (StringProperty, UniqueIdProperty, DateTimeProperty, RelationshipFrom)

if TYPE_CHECKING:
    from protrend.models.node import Node


class VersionNode:
    uid = UniqueIdProperty()
    name = StringProperty(required=True, unique_index=True)
    created = DateTimeProperty(default_now=True)

    __registered_relationships = {}

    @property
    def registered_relationships(self) -> Dict[str, RelationshipFrom]:
        return self.__registered_relationships

    @classmethod
    def register(cls, node: Type['Node']):

        setattr(cls, node.name(), RelationshipFrom(node, 'VERSIONING'))
        cls.__registered_relationships[node.name()] = getattr(cls, node.name())

    def get_children(self, node: Union[str, Type['Node']] = None):

        if not node:
            return {name: relationship.all()
                    for name, relationship in self.registered_relationships.keys()}

        if hasattr(node, 'name'):
            node = node.name()

        relationship = self.registered_relationships[node]
        return relationship.all()
