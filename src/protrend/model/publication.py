from neomodel import StringProperty, IntegerProperty, RelationshipTo

from protrend.utils.processors import rstrip, lstrip, lower_case, to_nan, to_int_str
from .base import BaseNode
from .relationships import REL_TYPE, BaseRelationship
from .utils import help_text


class Publication(BaseNode):
    entity = 'PUB'
    node_factors = {'pmid': [to_int_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    pmid = IntegerProperty(required=True, help_text=help_text.pmid)
    doi = StringProperty(max_length=250, help_text=help_text.doi)
    title = StringProperty(max_length=500, help_text=help_text.title)
    author = StringProperty(max_length=250, help_text=help_text.author)
    year = IntegerProperty(help_text=help_text.year)

    # relationships
    regulatory_family = RelationshipTo('.regulatory_family.RegulatoryFamily', REL_TYPE, model=BaseRelationship)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE, model=BaseRelationship)
    operon = RelationshipTo('.operon.Operon', REL_TYPE, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', REL_TYPE, model=BaseRelationship)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE, model=BaseRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE,
                                            model=BaseRelationship)
