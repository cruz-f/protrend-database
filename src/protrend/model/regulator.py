from neomodel import StringProperty, RelationshipTo, One

from .base import BaseNode, GeneMixIn, SequenceMixIn, PositionMixIn
from .relationships import REL_TYPE, SourceRelationship
from .utils import choices, help_text


class Regulator(BaseNode, GeneMixIn, SequenceMixIn, PositionMixIn):
    # base
    entity = 'REG'

    # properties
    mechanism = StringProperty(required=True, choices=choices.mechanism,
                               help_text=help_text.mechanism)

    # properties inherited from GeneMixIn, SequenceMixIn, PositionMixIn

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', REL_TYPE)
    publication = RelationshipTo('.publication.Publication', REL_TYPE)
    pathway = RelationshipTo('.pathway.Pathway', REL_TYPE)
    effector = RelationshipTo('.effector.Effector', REL_TYPE)
    regulatory_family = RelationshipTo('.regulatory_family.RegulatoryFamily', REL_TYPE, cardinality=One)
    organism = RelationshipTo('.organism.Organism', REL_TYPE, cardinality=One)
    gene = RelationshipTo('.gene.Gene', REL_TYPE)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
