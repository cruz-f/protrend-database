from neomodel import RelationshipTo, One, ZeroOrOne

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode, GeneMixIn, SequenceMixIn, PositionMixIn
from .relationships import REL_TYPE, SourceRelationship


class Gene(BaseNode, GeneMixIn, SequenceMixIn, PositionMixIn):
    # base
    entity = 'GEN'
    node_factors = {'uniprot_accession': [to_str, lower_case, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties inherited from GeneMixIn, SequenceMixIn, PositionMixIn

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', REL_TYPE)
    publication = RelationshipTo('.publication.Publication', REL_TYPE)
    pathway = RelationshipTo('.pathway.Pathway', REL_TYPE)
    operon = RelationshipTo('.operon.Operon', REL_TYPE)
    organism = RelationshipTo('.organism.Organism', REL_TYPE, cardinality=ZeroOrOne)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
