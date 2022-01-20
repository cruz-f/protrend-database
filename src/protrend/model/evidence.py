from neomodel import StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, RequiredNameMixIn
from .relationships import REL_TYPE, BaseRelationship
from .utils import help_text


class Evidence(BaseNode, RequiredNameMixIn):
    entity = 'EVI'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    description = StringProperty(help_text=help_text.generic_description)

    # relationships
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE, model=BaseRelationship)
    operon = RelationshipTo('.operon.Operon', REL_TYPE, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', REL_TYPE, model=BaseRelationship)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE, model=BaseRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE,
                                            model=BaseRelationship)
