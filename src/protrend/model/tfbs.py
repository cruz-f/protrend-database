from neomodel import StringProperty, RelationshipTo, IntegerProperty, One

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan
from .base import BaseNode, SequenceMixIn, PositionMixIn
from .relationships import REL_TYPE, SourceRelationship
from .utils import help_text


class TFBS(BaseNode, SequenceMixIn, PositionMixIn):
    entity = 'TBS'
    node_factors = {'site_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    site_hash = StringProperty(required=True, max_length=600)
    organism = StringProperty(required=True, max_length=100, help_text=help_text.organism_id)
    length = IntegerProperty(required=True, help_text=help_text.length)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', REL_TYPE)
    publication = RelationshipTo('.publication.Publication', REL_TYPE)
    data_organism = RelationshipTo('.organism.Organism', REL_TYPE, cardinality=One)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE)
    gene = RelationshipTo('.gene.Gene', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
