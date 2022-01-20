from neomodel import ArrayProperty, StringProperty, RelationshipTo

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, RequiredNameMixIn
from .relationships import REL_TYPE, SourceRelationship, BaseRelationship
from .utils import help_text


class Pathway(BaseNode, RequiredNameMixIn):
    entity = 'PTH'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    kegg_pathways = ArrayProperty(StringProperty(), help_text=help_text.kegg_pathways)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', REL_TYPE, model=BaseRelationship)
