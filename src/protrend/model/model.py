from neomodel import (StringProperty,
                      ArrayProperty,
                      RelationshipTo,
                      StructuredRel,
                      DateTimeProperty,
                      IntegerProperty,
                      One)

from protrend.model.node import Node
from protrend.utils.processors import to_str, rstrip, lstrip, lower_case, to_int_str, to_nan

# the main relationship type
REL_TYPE = 'HAS'


class SourceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()


class Source(Node):
    entity = 'SRC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    type = StringProperty(required=True)
    url = StringProperty()
    doi = StringProperty()
    authors = ArrayProperty(StringProperty())
    description = StringProperty()

    # relationships
    organism = RelationshipTo('Organism', REL_TYPE, model=SourceRelationship)
    pathway = RelationshipTo('Pathway', REL_TYPE, model=SourceRelationship)
    regulatory_family = RelationshipTo('RegulatoryFamily', REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('Regulator', REL_TYPE, model=SourceRelationship)
    operon = RelationshipTo('Operon', REL_TYPE, model=SourceRelationship)
    gene = RelationshipTo('Gene', REL_TYPE, model=SourceRelationship)
    tfbs = RelationshipTo('TFBS', REL_TYPE, model=SourceRelationship)
    effector = RelationshipTo('Effector', REL_TYPE, model=SourceRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE, model=SourceRelationship)


class Evidence(Node):
    entity = 'EVI'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    description = StringProperty()

    # relationships
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Publication(Node):
    entity = 'PUB'
    node_factors = {'pmid': [to_int_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    pmid = StringProperty(required=True)
    doi = StringProperty()
    title = StringProperty()
    author = StringProperty()
    year = StringProperty()

    # relationships
    regulatory_family = RelationshipTo('RegulatoryFamily', REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Pathway(Node):
    entity = 'PTH'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    kegg_pathways = ArrayProperty(StringProperty())

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)


class Effector(Node):
    entity = 'EFC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    kegg_compounds = ArrayProperty(StringProperty())

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class RegulatoryFamily(Node):
    entity = 'RFAM'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty(required=True)
    rfam = StringProperty()
    description = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)


class Operon(Node):
    entity = 'OPN'
    node_factors = {'operon_db_id': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    operon_db_id = StringProperty(required=True)
    name = StringProperty()
    function = StringProperty()
    genes = ArrayProperty(StringProperty(), required=True)
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, cardinality=One, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE, cardinality=One)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo('Organism', REL_TYPE, cardinality=One)
    gene = RelationshipTo('Gene', REL_TYPE)


class Organism(Node):
    entity = 'ORG'
    node_factors = {'ncbi_taxonomy': [to_int_str, lower_case, rstrip, lstrip, to_nan],
                    'name': [to_str, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    ncbi_taxonomy = IntegerProperty()
    species = StringProperty()
    strain = StringProperty()
    refseq_accession = StringProperty()
    refseq_ftp = StringProperty()
    genbank_accession = StringProperty()
    genbank_ftp = StringProperty()
    ncbi_assembly = IntegerProperty()
    assembly_accession = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    operon = RelationshipTo(Operon, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Regulator(Node):
    entity = 'REG'
    node_factors = {'uniprot_accession': [to_str, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty(required=True)
    uniprot_accession = StringProperty()
    name = StringProperty()
    synonyms = ArrayProperty(StringProperty())
    mechanism = StringProperty(required=True)
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = IntegerProperty()
    ncbi_protein = IntegerProperty()
    genbank_accession = StringProperty()
    refseq_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    pathway = RelationshipTo(Pathway, REL_TYPE)
    effector = RelationshipTo(Effector, REL_TYPE)
    regulatory_family = RelationshipTo(RegulatoryFamily, REL_TYPE, cardinality=One)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Gene(Node):
    entity = 'GEN'
    node_factors = {'uniprot_accession': [to_str, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty(required=True)
    uniprot_accession = StringProperty()
    name = StringProperty()
    synonyms = ArrayProperty(StringProperty())
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = IntegerProperty()
    ncbi_protein = IntegerProperty()
    genbank_accession = StringProperty()
    refseq_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    pathway = RelationshipTo(Pathway, REL_TYPE)
    operon = RelationshipTo(Operon, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class TFBS(Node):
    entity = 'TBS'
    node_factors = {'site_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    site_hash = StringProperty(required=True)
    organism = StringProperty(required=True)
    sequence = StringProperty(required=True)
    strand = StringProperty(required=True)
    start = IntegerProperty(required=True)
    stop = IntegerProperty(required=True)
    length = IntegerProperty(required=True)

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    data_organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE)
    gene = RelationshipTo(Gene, REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class RegulatoryInteraction(Node):
    entity = 'RIN'
    node_factors = {'regulatory_interaction_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    regulatory_interaction_hash = StringProperty(required=True)
    organism = StringProperty(required=True)
    regulator = StringProperty(required=True)
    gene = ArrayProperty(required=True)
    tfbs = ArrayProperty()
    effector = StringProperty()
    regulatory_effect = StringProperty(required=True)

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    data_effector = RelationshipTo(Effector, REL_TYPE, cardinality=One)
    data_organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    data_regulator = RelationshipTo(Regulator, REL_TYPE, cardinality=One)
    data_gene = RelationshipTo(Gene, REL_TYPE, cardinality=One)
    data_tfbs = RelationshipTo(TFBS, REL_TYPE, cardinality=One)
