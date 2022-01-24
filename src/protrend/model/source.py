from neomodel import StringProperty, RelationshipTo, ArrayProperty

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, NameMixIn
from .relationships import SourceRelationship, SOURCE_REL_TYPE
from .utils import help_text, choices


class Source(BaseNode, NameMixIn):
    entity = 'SRC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    type = StringProperty(required=True, choices=choices.data_source_type,
                          help_text=help_text.data_source_type)
    url = StringProperty(max_length=300, help_text=help_text.url)
    doi = StringProperty(max_length=250, help_text=help_text.doi)
    authors = ArrayProperty(StringProperty(), help_text=help_text.source_author)
    description = StringProperty(help_text=help_text.generic_description)

    # relationships
    organism = RelationshipTo('.organism.Organism', SOURCE_REL_TYPE, model=SourceRelationship)
    pathway = RelationshipTo('.pathway.Pathway', SOURCE_REL_TYPE, model=SourceRelationship)
    regulatory_family = RelationshipTo('.regulatory_family.RegulatoryFamily', SOURCE_REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('.regulator.Regulator', SOURCE_REL_TYPE, model=SourceRelationship)
    operon = RelationshipTo('.operon.Operon', SOURCE_REL_TYPE, model=SourceRelationship)
    gene = RelationshipTo('.gene.Gene', SOURCE_REL_TYPE, model=SourceRelationship)
    tfbs = RelationshipTo('.tfbs.TFBS', SOURCE_REL_TYPE, model=SourceRelationship)
    effector = RelationshipTo('.effector.Effector', SOURCE_REL_TYPE, model=SourceRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', SOURCE_REL_TYPE,
                                            model=SourceRelationship)
