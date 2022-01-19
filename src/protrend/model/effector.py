from neomodel import ArrayProperty, StringProperty, RelationshipTo

from .base import BaseNode, RequiredNameMixIn
from .relationships import REL_TYPE, SourceRelationship
from .utils import help_text


class Effector(BaseNode, RequiredNameMixIn):
    entity = 'EFC'

    # properties
    kegg_compounds = ArrayProperty(StringProperty(), help_text=help_text.kegg_compounds)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
