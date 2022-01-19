from neomodel import StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, RequiredNameMixIn
from .relationships import REL_TYPE
from .utils import help_text


class Evidence(BaseNode, RequiredNameMixIn):
    entity = 'EVI'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    description = StringProperty(help_text=help_text.generic_description)

    # relationships
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE)
    operon = RelationshipTo('.operon.Operon', REL_TYPE)
    gene = RelationshipTo('.gene.Gene', REL_TYPE)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
