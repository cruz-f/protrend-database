from neomodel import (StringProperty,
                      ArrayProperty,
                      RelationshipTo,
                      StructuredRel,
                      DateTimeProperty,
                      UniqueIdProperty,
                      IntegerProperty,
                      One,
                      OneOrMore)

from protrend.model.node import Node

# the main relationship type
REL_TYPE = 'HAS'


class Source(Node):
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
    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()


class Evidence(Node):
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
    # properties
    doi = StringProperty(required=True)
    pmid = StringProperty(required=True)
    title = StringProperty(required=True)
    author = StringProperty(required=True)
    year = StringProperty(required=True)

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
    # properties
    name = StringProperty(required=True)
    species = StringProperty(required=True)
    strain = StringProperty()
    family = StringProperty()
    phylum = StringProperty()
    ncbi_taxonomy = StringProperty(required=True)
    refseq_accession = StringProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
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
    # properties
    name = StringProperty(required=True)
    kegg_pathway = StringProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)
    gene = RelationshipTo('Gene', REL_TYPE)


class RegulatoryFamily(Node):
    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty(required=True)
    rfam = StringProperty()
    description = StringProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    regulator = RelationshipTo('Regulator', REL_TYPE)


class OperonRelationship:
    # base
    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    operon = StringProperty()


class Regulator(Node):
    # properties
    locus_tag = StringProperty(required=True)
    name = StringProperty(required=True)
    rfam = StringProperty()
    synonyms = ArrayProperty(StringProperty())
    mechanism = StringProperty(required=True)
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = StringProperty()
    ncbi_protein = StringProperty()
    ncbi_accession = StringProperty()
    uniprot_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    position_left = IntegerProperty()
    position_right = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
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
    # properties
    name = StringProperty(required=True)
    promoters = ArrayProperty(StringProperty())
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    strand = StringProperty()
    first_gene_position_left = IntegerProperty()
    first_gene_position_right = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
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
    # properties
    name = StringProperty(required=True)
    sequence = StringProperty()
    strand = StringProperty()
    position_left = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo('Gene', REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)


class Gene(Node):
    # properties
    locus_tag = StringProperty(required=True)
    name = StringProperty(required=True)
    synonyms = ArrayProperty(StringProperty())
    function = StringProperty()
    description = StringProperty()
    ncbi_gene = StringProperty()
    ncbi_protein = StringProperty()
    refseq_accession = StringProperty()
    uniprot_accession = StringProperty()
    sequence = StringProperty()
    strand = StringProperty()
    position_left = IntegerProperty()
    position_right = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
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
    # properties
    sequence = StringProperty(required=True)
    strand = StringProperty()
    position_left = IntegerProperty()
    position_right = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE, model=OperonRelationship)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    regulon = RelationshipTo('Regulon', REL_TYPE, model=OperonRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE, model=OperonRelationship)


class Effector(Node):
    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty()
    kegg_metabolite = StringProperty()

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE)
    regulator = RelationshipTo(Regulator, REL_TYPE, cardinality=OneOrMore)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', REL_TYPE)


class Regulon(Node):
    # properties
    regulators = ArrayProperty(StringProperty(), required=True)
    operons = ArrayProperty(StringProperty(), required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    type = StringProperty(required=True)

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator = RelationshipTo(Regulator, REL_TYPE, cardinality=OneOrMore)
    operon = RelationshipTo(Operon, REL_TYPE, cardinality=OneOrMore)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, REL_TYPE, model=OperonRelationship)


class RegulatoryInteraction(Node):
    # properties
    regulator = StringProperty(required=True)
    operon = StringProperty(required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    effectors = ArrayProperty(StringProperty())
    regulatory_effect = StringProperty(required=True)

    # relationships
    source = RelationshipTo(Source, REL_TYPE, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, REL_TYPE)
    publication = RelationshipTo(Publication, REL_TYPE)
    organism = RelationshipTo(Organism, REL_TYPE, cardinality=One)
    regulator_in = RelationshipTo(Regulator, REL_TYPE, cardinality=One)
    operon_out = RelationshipTo(Operon, REL_TYPE, cardinality=One)
    gene = RelationshipTo(Gene, REL_TYPE, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, REL_TYPE, model=OperonRelationship)
    effector = RelationshipTo(Effector, REL_TYPE)
