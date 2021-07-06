from neomodel import (StringProperty)

from protrend.model.node import Node


class Version(Node):
    property_as_id = 'name'

    name = StringProperty(required=True, unique_index=True)

    __abstract_node__ = True
