from neomodel import ArrayProperty, StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, NameMixIn
from .relationships import SOURCE_REL_TYPE, BASE_REL_TYPE, SourceRelationship, BaseRelationship
from .utils import help_text


class Effector(BaseNode, NameMixIn):
    entity = 'EFC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    kegg_compounds = ArrayProperty(StringProperty(), help_text=help_text.kegg_compounds)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('.regulator.Regulator', BASE_REL_TYPE, model=BaseRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', BASE_REL_TYPE,
                                            model=BaseRelationship)
