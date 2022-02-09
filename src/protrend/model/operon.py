from neomodel import StringProperty, ArrayProperty, RelationshipTo, One, IntegerProperty

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan
from .base import BaseNode
from .relationships import SourceRelationship, BaseRelationship, SOURCE_REL_TYPE, BASE_REL_TYPE
from .utils import help_text, choices


class Operon(BaseNode):
    # base
    entity = 'OPN'
    node_factors = {'operon_db_id': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    operon_db_id = StringProperty(required=True, unique_index=True, max_length=50, help_text=help_text.operon_db_id)
    operon_db_id_factor = StringProperty(required=True, unique_index=True, max_length=50,
                                         help_text=help_text.operon_db_id)
    name = StringProperty(max_length=50, help_text=help_text.operon_name)
    function = StringProperty(max_length=250, help_text=help_text.operon_function)
    genes = ArrayProperty(StringProperty(), required=True, help_text=help_text.operon_genes)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, cardinality=One, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', BASE_REL_TYPE, cardinality=One, model=BaseRelationship)
    publication = RelationshipTo('.publication.Publication', BASE_REL_TYPE, model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, cardinality=One, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, model=BaseRelationship)
