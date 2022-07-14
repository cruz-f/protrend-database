from neomodel import StringProperty, IntegerProperty, RelationshipTo, ZeroOrOne

from protrend.utils.processors import rstrip, lstrip, lower_case, to_nan, to_str
from .base import BaseNode
from .relationships import BaseRelationship, BASE_REL_TYPE, SOURCE_REL_TYPE, SourceRelationship
from .utils import help_text, choices


class PromoterRegion(BaseNode):
    entity = 'PRO'
    node_factors = {'locus_tag': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty(required=True, unique_index=True, max_length=100, help_text=help_text.locus_tag)
    locus_tag_factor = StringProperty(required=True, unique_index=True, max_length=100,
                                      help_text=help_text.locus_tag)
    gene = StringProperty(required=True, max_length=100, help_text=help_text.gene_id)
    sequence = StringProperty(help_text=help_text.promoter_sequence)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, cardinality=ZeroOrOne, model=SourceRelationship)
    data_gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)
