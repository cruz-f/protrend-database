from neomodel import StringProperty, RelationshipTo, IntegerProperty, ZeroOrOne

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan
from .base import BaseNode
from .relationships import SourceRelationship, BaseRelationship, BASE_REL_TYPE, SOURCE_REL_TYPE
from .utils import help_text, choices


class TFBS(BaseNode):
    entity = 'TBS'
    node_factors = {'site_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    site_hash = StringProperty(required=True, unique_index=True, max_length=600)
    site_hash_factor = StringProperty(required=True, unique_index=True, max_length=600)
    organism = StringProperty(required=True, max_length=100, help_text=help_text.organism_id)
    sequence = StringProperty(help_text=help_text.tfbs_sequence)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)
    length = IntegerProperty(required=True, help_text=help_text.length)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', BASE_REL_TYPE, model=BaseRelationship)
    publication = RelationshipTo('.publication.Publication', BASE_REL_TYPE, model=BaseRelationship)
    data_organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)
    regulator = RelationshipTo('.regulator.Regulator', BASE_REL_TYPE, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, model=BaseRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', BASE_REL_TYPE,
                                            model=BaseRelationship)
    motif = RelationshipTo('.motif.Motif', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)