from neomodel import StringProperty, RelationshipTo, One

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan
from .base import BaseNode
from .relationships import REL_TYPE, SourceRelationship
from .utils import help_text, choices


class RegulatoryInteraction(BaseNode):
    entity = 'RIN'
    node_factors = {'regulatory_interaction_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    regulatory_interaction_hash = StringProperty(required=True, max_length=600)
    organism = StringProperty(required=True, max_length=100, help_text=help_text.organism_id)
    regulator = StringProperty(required=True, max_length=100, help_text=help_text.regulator_id)
    gene = StringProperty(required=True, max_length=100, help_text=help_text.gene_id)
    tfbs = StringProperty(max_length=100, help_text=help_text.tfbs_id)
    effector = StringProperty(max_length=100, help_text=help_text.effector_id)
    regulatory_effect = StringProperty(required=True, choices=choices.regulatory_effect,
                                       help_text=help_text.regulatory_effect)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', REL_TYPE)
    publication = RelationshipTo('.publication.Publication', REL_TYPE)
    data_effector = RelationshipTo('.effector.Effector', REL_TYPE, cardinality=One)
    data_organism = RelationshipTo('.organism.Organism', REL_TYPE, cardinality=One)
    data_regulator = RelationshipTo('.regulator.Regulator', REL_TYPE, cardinality=One)
    data_gene = RelationshipTo('.gene.Gene', REL_TYPE, cardinality=One)
    data_tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE, cardinality=One)
