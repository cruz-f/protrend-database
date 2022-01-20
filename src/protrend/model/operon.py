from neomodel import StringProperty, ArrayProperty, RelationshipTo, One

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan
from .base import BaseNode, PositionMixIn
from .relationships import REL_TYPE, SourceRelationship, BaseRelationship
from .utils import help_text


class Operon(BaseNode, PositionMixIn):
    # base
    entity = 'OPN'
    node_factors = {'operon_db_id': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties inherited from PositionMixIn

    # properties
    operon_db_id = StringProperty(required=True, max_length=50, help_text=help_text.operon_db_id)
    name = StringProperty(max_length=50, help_text=help_text.operon_name)
    function = StringProperty(max_length=250, help_text=help_text.operon_function)
    genes = ArrayProperty(StringProperty(), required=True, help_text=help_text.operon_genes)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, cardinality=One, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', REL_TYPE, cardinality=One, model=BaseRelationship)
    publication = RelationshipTo('.publication.Publication', REL_TYPE, model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', REL_TYPE, cardinality=One, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', REL_TYPE, model=BaseRelationship)
