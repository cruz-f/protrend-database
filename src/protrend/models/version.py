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

        setattr(cls, node.cls_name(), RelationshipFrom(node, 'VERSIONING'))
        cls.__registered_relationships[node.cls_name()] = getattr(cls, node.cls_name())

    def get_children(self,
                     node_type: str = None,
                     to: str = 'list') -> Union[list, dict]:

        if node_type:

            if to == 'dict':

                return {node_type: self.registered_relationships[node_type].all()}

            elif to == 'list':

                return self.registered_relationships[node_type].all()

            else:
                raise ValueError(f'{to} output format not supported')

        else:

            if to == 'dict':
                return {name: relationship.all()
                        for name, relationship in self.registered_relationships.items()}

            elif to == 'list':
                nodes = []
                for relationship in self.registered_relationships.values():
                    nodes.extend(relationship.all())

                return nodes

            else:
                raise ValueError(f'{to} output format not supported')
