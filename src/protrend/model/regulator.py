from neomodel import StringProperty, RelationshipTo, ZeroOrOne, ArrayProperty, IntegerProperty

from protrend.utils.processors import to_str, lower_case, rstrip, lstrip, to_nan
from .base import BaseNode
from .relationships import SourceRelationship, BaseRelationship, BASE_REL_TYPE, SOURCE_REL_TYPE
from .utils import choices, help_text


class Regulator(BaseNode):
    # base
    entity = 'REG'
    node_factors = {'uniprot_accession': [to_str, lower_case, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty(required=True, unique_index=True, max_length=100, help_text=help_text.locus_tag)
    locus_tag_factor = StringProperty(required=True, unique_index=True, max_length=100,
                                      help_text=help_text.required_name)
    uniprot_accession = StringProperty(unique_index=True, max_length=50, help_text=help_text.uniprot_accession)
    uniprot_accession_factor = StringProperty(unique_index=True, max_length=50, help_text=help_text.uniprot_accession)
    name = StringProperty(max_length=50, help_text=help_text.gene_name)
    synonyms = ArrayProperty(StringProperty(), help_text=help_text.synonyms)
    function = StringProperty(help_text=help_text.function)
    description = StringProperty(help_text=help_text.description)
    mechanism = StringProperty(required=True, choices=choices.mechanism,
                               help_text=help_text.mechanism)
    ncbi_gene = IntegerProperty(max_length=50, help_text=help_text.ncbi_gene)
    ncbi_protein = IntegerProperty(max_length=50, help_text=help_text.ncbi_protein)
    genbank_accession = StringProperty(max_length=50, help_text=help_text.genbank_accession)
    refseq_accession = StringProperty(max_length=50, help_text=help_text.refseq_accession)
    sequence = StringProperty(help_text=help_text.sequence)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)

    # relationships
    data_source = RelationshipTo('.source.Source', SOURCE_REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo('.evidence.Evidence', BASE_REL_TYPE, model=BaseRelationship)
    publication = RelationshipTo('.publication.Publication', BASE_REL_TYPE, model=BaseRelationship)
    pathway = RelationshipTo('.pathway.Pathway', BASE_REL_TYPE, model=BaseRelationship)
    effector = RelationshipTo('.effector.Effector', BASE_REL_TYPE, model=BaseRelationship)
    regulatory_family = RelationshipTo('.regulatory_family.RegulatoryFamily', BASE_REL_TYPE, cardinality=ZeroOrOne,
                                       model=BaseRelationship)
    organism = RelationshipTo('.organism.Organism', BASE_REL_TYPE, cardinality=ZeroOrOne, model=BaseRelationship)
    gene = RelationshipTo('.gene.Gene', BASE_REL_TYPE, model=BaseRelationship)
    tfbs = RelationshipTo('.tfbs.TFBS', BASE_REL_TYPE, model=BaseRelationship)
    regulatory_interaction = RelationshipTo('.regulatory_interaction.RegulatoryInteraction', BASE_REL_TYPE,
                                            model=BaseRelationship)
