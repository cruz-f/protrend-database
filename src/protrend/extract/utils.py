from typing import List, Type, Union

from neomodel import StructuredNode

from protrend.models.node import Node
from protrend.models.version import VersionNode


def parsing_spider_arguments(argument: str):

    if not argument:
        return

    children = argument.split(',')

    return tuple(children)


class NodeRelationshipMap:

    def __init__(self,
                 node: Node,
                 relationship_name: str,
                 to_node_cls: Union[Type[StructuredNode, Type[Node]]],
                 to_node_attr: str,
                 to: List[str, int]):

        """

        NodeRelationshipMap holds the link between nodes and identifiers of the children node, as connections can
        only be made after the creation of all nodes. Note that, different node types are created at different moments

        :param node:
        :param relationship_name:
        :param to_node_cls:
        :param to_node_attr:
        :param to:
        """

        self.node = node
        self.relationship_name = relationship_name
        self.to_node_cls = to_node_cls
        self.to_node_attr = to_node_attr
        self.to = to

    def connect(self, version: VersionNode):

        for node_id in self.to:
            to_node = self.to_node_cls.get_by_version(attr=self.to_node_attr, value=node_id, version=version)

            if to_node:
                rel = getattr(self.node, self.relationship_name)
                rel.connect(to_node)
