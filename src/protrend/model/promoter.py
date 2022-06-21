from neomodel import StringProperty, IntegerProperty, RelationshipTo

from protrend.utils.processors import rstrip, lstrip, lower_case, to_nan, to_int_str
from .base import BaseNode
from .relationships import BaseRelationship, BASE_REL_TYPE
from .utils import help_text


class Promoter(BaseNode):
    entity = 'PRO'
    node_factors = {'pmid': [to_int_str, lower_case, rstrip, lstrip, to_nan]} # compor

    # properties
    promoter_sequence = StringProperty(help_text=help_text.promoter_sequence)
    locus_tag = StringProperty(required=True, unique_index=True, max_length=100, help_text=help_text.locus_tag)
    uniprot_accession = StringProperty(unique_index=True, max_length=50, help_text=help_text.uniprot_accession)
    start = IntegerProperty(help_text=help_text.start) # uso o start e stop já existente ou crio específicos para promoters?
    end = IntegerProperty(help_text=help_text.stop)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, model=BaseRelationship)
