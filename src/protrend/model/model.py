from neomodel import (StringProperty,
                      ArrayProperty,
                      RelationshipTo,
                      StructuredRel,
                      DateTimeProperty,
                      IntegerProperty,
                      One,
                      OneOrMore)

from protrend.model.node import Node
from protrend.transform.processors import to_str, rstrip, lstrip, lower_case, to_int_str, to_nan

# the main relationship type
REL_TYPE = 'HAS'


class Source(Node):
    entity = 'SRC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    type = StringProperty(required=True)
    url = StringProperty()
    doi = StringProperty()
    authors = ArrayProperty(StringProperty(), required=True)
    description = StringProperty()

    # relationships
    organism = RelationshipTo('Organism', REL_TYPE)
    pathway = RelationshipTo('Pathway', REL_TYPE)
    regulatory_family = RelationshipTo('RegulatoryFamily', REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    promoter = RelationshipTo('Promoter', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    effector = RelationshipTo('Effector', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class SourceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()


class Evidence(Node):
    entity = 'EVI'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    description = StringProperty()

    # relationships
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    promoter = RelationshipTo('Promoter', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Publication(Node):
    entity = 'PUB'
    node_factors = {'pmid': [to_int_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    doi = StringProperty()
    pmid = StringProperty()
    title = StringProperty()
    author = StringProperty()
    year = StringProperty()

    # relationships
    organism = RelationshipTo('Organism', REL_TYPE)
    pathway = RelationshipTo('Pathway', REL_TYPE)
    regulatory_family = RelationshipTo('RegulatoryFamily', REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    promoter = RelationshipTo('Promoter', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    effector = RelationshipTo('Effector', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Organism(Node):
    entity = 'ORG'
    node_factors = {'ncbi_taxonomy': [to_int_str, lower_case, rstrip, lstrip, to_nan],
                    'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    species = StringProperty()
    strain = StringProperty()
    ncbi_taxonomy = IntegerProperty()
    refseq_accession = StringProperty()
    refseq_ftp = StringProperty()
    genbank_accession = StringProperty()
    genbank_ftp = StringProperty()
    ncbi_assembly = IntegerProperty()
    assembly_accession = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    promoter = RelationshipTo('Promoter', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    effector = RelationshipTo('Effector', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Pathway(Node):
    entity = 'PTH'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    kegg_pathways = ArrayProperty(StringProperty())

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)


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


class OperonRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    operon = StringProperty()


class Regulator(Node):
    entity = 'REG'
    node_factors = {'uniprot_accession': [to_str, lower_case, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty()
    name = StringProperty(required=True)
    rfam = StringProperty()
    synonyms = ArrayProperty(StringProperty())
    mechanism = StringProperty(required=True)
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = IntegerProperty()
    ncbi_protein = IntegerProperty()
    genbank_accession = StringProperty()
    refseq_accession = StringProperty()
    uniprot_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    pathway = RelationshipTo(Pathway, REL_TYPE)
    regulatory_family = RelationshipTo(RegulatoryFamily, REL_TYPE)
    operon = RelationshipTo('Operon', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE, model=OperonRelationship)
    tfbs = RelationshipTo('TFBS', REL_TYPE, model=OperonRelationship)
    effector = RelationshipTo('Effector', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Operon(Node):
    entity = 'OPN'
    node_factors = {'operon_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty()
    promoters = ArrayProperty(StringProperty())
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()
    operon_hash = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE)
    promoter = RelationshipTo('Promoter', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE, cardinality=OneOrMore)
    tfbs = RelationshipTo('TFBS', REL_TYPE)
    regulon = RelationshipTo('Regulon', REL_TYPE)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Promoter(Node):
    entity = 'PRO'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    sequence = StringProperty()
    strand = StringProperty()
    start = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo('Gene', REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)


class Gene(Node):
    entity = 'GEN'
    node_factors = {'uniprot_accession': [to_str, lower_case, rstrip, lstrip, to_nan],
                    'locus_tag': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    locus_tag = StringProperty()
    name = StringProperty(required=True)
    synonyms = ArrayProperty(StringProperty())
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = IntegerProperty()
    ncbi_protein = IntegerProperty()
    genbank_accession = StringProperty()
    refseq_accession = StringProperty()
    uniprot_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    pathway = RelationshipTo(Pathway, REL_TYPE)
    regulator = RelationshipTo(Regulator, REL_TYPE, model=OperonRelationship)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    promoter = RelationshipTo(Promoter, REL_TYPE, model=OperonRelationship)
    tfbs = RelationshipTo('TFBS', REL_TYPE, model=OperonRelationship)
    regulon = RelationshipTo('Regulon', REL_TYPE, model=OperonRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE, model=OperonRelationship)


class TFBS(Node):
    entity = 'TBS'
    node_factors = {'site_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    sequence = StringProperty(required=True)
    strand = StringProperty()
    start = IntegerProperty()
    stop = IntegerProperty()
    length = IntegerProperty()
    site_hash = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE, model=OperonRelationship)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    regulon = RelationshipTo('Regulon', REL_TYPE, model=OperonRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE, model=OperonRelationship)


class Effector(Node):
    entity = 'EFC'
    node_factors = {'name': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty()
    kegg_compounds = ArrayProperty(StringProperty())

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE)
    regulator = RelationshipTo(Regulator, REL_TYPE, cardinality=OneOrMore)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Regulon(Node):
    entity = 'REN'

    # properties
    regulators = ArrayProperty(StringProperty(), required=True)
    operons = ArrayProperty(StringProperty(), required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    type = StringProperty(required=True)

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE, cardinality=OneOrMore)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=OneOrMore)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, REL_TYPE, model=OperonRelationship)


class RegulatoryInteraction(Node):
    entity = 'RIN'
    node_factors = {'regulatory_interaction_hash': [to_str, lower_case, rstrip, lstrip, to_nan]}

    # properties
    regulator = StringProperty(required=True)
    operon = StringProperty(required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    effectors = ArrayProperty(StringProperty())
    regulatory_effect = StringProperty()
    regulatory_interaction_hash = StringProperty()

    # relationships
    data_source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator_in = RelationshipTo(Regulator, REL_TYPE, cardinality=One)
    operon_out = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, REL_TYPE, model=OperonRelationship)
    effector = RelationshipTo(Effector, REL_TYPE)
