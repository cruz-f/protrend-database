from neomodel import StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode
from .relationships import SourceRelationship, BaseRelationship, BASE_REL_TYPE, SOURCE_REL_TYPE
from .utils import help_text, choices


class RegulatoryFamily(BaseNode):
    entity = 'RFAM'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)
    name_factor = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)
    mechanism = StringProperty(required=True, choices=choices.mechanism,
                               help_text=help_text.mechanism)
    rfam = StringProperty(max_length=100, help_text=help_text.rfam)
    description = StringProperty(help_text=help_text.generic_description)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo('.publication.Publication', BASE_REL_TYPE, model=BaseRelationship)
    regulator = RelationshipTo('.regulator.Regulator', BASE_REL_TYPE, model=BaseRelationship)
