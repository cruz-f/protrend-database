from typing import Union, Dict

from neomodel import (StringProperty, UniqueIdProperty, DateTimeProperty, StructuredNode,
                      RelationshipTo)


class Version(StructuredNode):

    uid = UniqueIdProperty()
    name = StringProperty(required=True, unique_index=True)
    created = DateTimeProperty(default_now=True)

    __abstract_node__ = True

    @property
    def versioned_nodes(self) -> Dict[str, RelationshipTo]:
        return {}

    def get_versioned_nodes(self,
                            node_type: str = None,
                            to: str = 'list') -> Union[list, dict]:

        if node_type:

            if to == 'dict':

                return {node_type: self.versioned_nodes[node_type].all()}

            elif to == 'list':

                return self.versioned_nodes[node_type].all()

            else:
                raise ValueError(f'{to} output format not supported')

        else:

            if to == 'dict':
                return {name: relationship.all()
                        for name, relationship in self.versioned_nodes.items()}

            elif to == 'list':
                nodes = []
                for relationship in self.versioned_nodes.values():
                    nodes.extend(relationship.all())

                return nodes

            else:
                raise ValueError(f'{to} output format not supported')
