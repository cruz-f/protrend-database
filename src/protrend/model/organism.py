from neomodel import StringProperty, IntegerProperty, RelationshipTo

from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_nan, to_int_str
from .base import BaseNode
from .relationships import REL_TYPE, SourceRelationship
from .utils import help_text


class Organism(BaseNode):
    # base
    entity = 'ORG'
    node_factors = {'ncbi_taxonomy': [to_int_str, lower_case, rstrip, lstrip, to_nan],
                    'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True, max_length=200, help_text=help_text.organism_name)
    ncbi_taxonomy = IntegerProperty(help_text=help_text.ncbi_taxonomy)
    species = StringProperty(max_length=150, help_text=help_text.species)
    strain = StringProperty(max_length=150, help_text=help_text.strain)
    refseq_accession = StringProperty(max_length=50, help_text=help_text.refseq_accession)
    refseq_ftp = StringProperty(max_length=250, help_text=help_text.refseq_ftp)
    genbank_accession = StringProperty(max_length=50, help_text=help_text.genbank_accession)
    genbank_ftp = StringProperty(max_length=250, help_text=help_text.genbank_ftp)
    ncbi_assembly = IntegerProperty(help_text=help_text.ncbi_assembly)
    assembly_accession = StringProperty(max_length=50, help_text=help_text.assembly_accession)

    # relationships
    data_source = RelationshipTo('.source.Source', REL_TYPE, model=SourceRelationship)
    operon = RelationshipTo('.operon.Operon', REL_TYPE)
    regulator = RelationshipTo('.regulator.Regulator', REL_TYPE)
    gene = RelationshipTo('.gene.Gene', REL_TYPE)
    tfbs = RelationshipTo('.tfbs.TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', REL_TYPE)
