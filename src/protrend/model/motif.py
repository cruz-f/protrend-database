from neomodel import StringProperty, IntegerProperty, RelationshipTo, ZeroOrOne

from protrend.utils.processors import rstrip, lstrip, lower_case, to_nan, to_int_str
from .base import BaseNode
from .relationships import BaseRelationship, BASE_REL_TYPE, SOURCE_REL_TYPE, SourceRelationship
from .utils import help_text, choices


class Motif(BaseNode):
    entity = 'TFM'
    node_factors = {'pmid': [to_int_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    motif_sequence = StringProperty(help_text=help_text.motif_sequence)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
    score = IntegerProperty(help_text=help_text.score)
    start = IntegerProperty(help_text=help_text.start)
    end = IntegerProperty(help_text=help_text.stop)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, cardinality=ZeroOrOne, model=SourceRelationship)
    tfbs = RelationshipTo('.tfbs.TFBS', BASE_REL_TYPE, model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)
    regulator = RelationshipTo('.regulator.Regulator', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)