from neomodel import ArrayProperty, StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode
from .relationships import SourceRelationship, BaseRelationship, SOURCE_REL_TYPE, BASE_REL_TYPE
from .utils import help_text


class Pathway(BaseNode):
    entity = 'PTH'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)
    name_factor = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)
    kegg_pathways = ArrayProperty(StringProperty(), help_text=help_text.kegg_pathways)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('.regulator.Regulator', BASE_REL_TYPE, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, model=BaseRelationship)
