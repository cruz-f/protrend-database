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
HAS = 'HAS'


class Source(Node):
    # properties
    name = StringProperty(required=True)
    type = StringProperty(required=True)
    url = StringProperty()
    doi = StringProperty()
    authors = ArrayProperty(StringProperty(), required=True)
    description = StringProperty()

    # relationships
    organism = RelationshipTo('Organism', HAS)
    pathway = RelationshipTo('Pathway', HAS)
    regulatory_family = RelationshipTo('RegulatoryFamily', HAS)
    regulator = RelationshipTo('Regulator', HAS)
    operon = RelationshipTo('Operon', HAS)
    promoter = RelationshipTo('Promoter', HAS)
    gene = RelationshipTo('Gene', HAS)
    tfbs = RelationshipTo('TFBS', HAS)
    effector = RelationshipTo('Effector', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


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
    regulator = RelationshipTo('Regulator', HAS)
    operon = RelationshipTo('Operon', HAS)
    promoter = RelationshipTo('Promoter', HAS)
    gene = RelationshipTo('Gene', HAS)
    tfbs = RelationshipTo('TFBS', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


class Publication(Node):
    # properties
    doi = StringProperty(required=True)
    pmid = StringProperty(required=True)
    title = StringProperty(required=True)
    author = StringProperty(required=True)
    year = StringProperty(required=True)

    # relationships
    organism = RelationshipTo('Organism', HAS)
    pathway = RelationshipTo('Pathway', HAS)
    regulatory_family = RelationshipTo('RegulatoryFamily', HAS)
    regulator = RelationshipTo('Regulator', HAS)
    operon = RelationshipTo('Operon', HAS)
    promoter = RelationshipTo('Promoter', HAS)
    gene = RelationshipTo('Gene', HAS)
    tfbs = RelationshipTo('TFBS', HAS)
    effector = RelationshipTo('Effector', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


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
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    publication = RelationshipTo(Publication, HAS)
    regulator = RelationshipTo('Regulator', HAS)
    operon = RelationshipTo('Operon', HAS)
    promoter = RelationshipTo('Promoter', HAS)
    gene = RelationshipTo('Gene', HAS)
    tfbs = RelationshipTo('TFBS', HAS)
    effector = RelationshipTo('Effector', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


class Pathway(Node):
    # properties
    name = StringProperty(required=True)
    kegg_pathway = StringProperty()

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    publication = RelationshipTo(Publication, HAS)
    regulator = RelationshipTo('Regulator', HAS)
    gene = RelationshipTo('Gene', HAS)


class RegulatoryFamily(Node):
    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty(required=True)
    rfam = StringProperty()
    description = StringProperty()

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    publication = RelationshipTo(Publication, HAS)
    regulator = RelationshipTo('Regulator', HAS)


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
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    pathway = RelationshipTo(Pathway, HAS)
    regulatory_family = RelationshipTo(RegulatoryFamily, HAS)
    operon = RelationshipTo('Operon', HAS)
    gene = RelationshipTo('Gene', HAS, model=OperonRelationship)
    tfbs = RelationshipTo('TFBS', HAS, model=OperonRelationship)
    effector = RelationshipTo('Effector', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


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
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    regulator = RelationshipTo(Regulator, HAS)
    promoter = RelationshipTo('Promoter', HAS)
    gene = RelationshipTo('Gene', HAS, cardinality=OneOrMore)
    tfbs = RelationshipTo('TFBS', HAS)
    regulon = RelationshipTo('Regulon', HAS)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


class Promoter(Node):
    # properties
    name = StringProperty(required=True)
    sequence = StringProperty()
    strand = StringProperty()
    position_left = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    operon = RelationshipTo(Operon, HAS, cardinality=One)
    gene = RelationshipTo('Gene', HAS, cardinality=OneOrMore, model=OperonRelationship)


class Gene(Node):
    # properties
    locus_tag = StringProperty(required=True)
    name = StringProperty(required=True)
    synonyms = ArrayProperty(StringProperty())
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
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    pathway = RelationshipTo(Pathway, HAS)
    regulator = RelationshipTo(Regulator, HAS, model=OperonRelationship)
    operon = RelationshipTo(Operon, HAS, cardinality=One)
    promoter = RelationshipTo(Promoter, HAS, model=OperonRelationship)
    tfbs = RelationshipTo('TFBS', HAS, model=OperonRelationship)
    regulon = RelationshipTo('Regulon', HAS, model=OperonRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS, model=OperonRelationship)


class TFBS(Node):
    # properties
    sequence = StringProperty(required=True)
    strand = StringProperty()
    position_left = IntegerProperty()
    position_right = IntegerProperty()

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    regulator = RelationshipTo(Regulator, HAS, model=OperonRelationship)
    operon = RelationshipTo(Operon, HAS, cardinality=One)
    gene = RelationshipTo(Gene, HAS, cardinality=OneOrMore, model=OperonRelationship)
    regulon = RelationshipTo('Regulon', HAS, model=OperonRelationship)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS, model=OperonRelationship)


class Effector(Node):
    # properties
    name = StringProperty(required=True)
    mechanism = StringProperty()
    kegg_metabolite = StringProperty()

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS)
    regulator = RelationshipTo(Regulator, HAS, cardinality=OneOrMore)
    regulatory_interaction = RelationshipTo('RegulatoryInteraction', HAS)


class Regulon(Node):
    # properties
    regulators = ArrayProperty(StringProperty(), required=True)
    operons = ArrayProperty(StringProperty(), required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    type = StringProperty(required=True)

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    regulator = RelationshipTo(Regulator, HAS, cardinality=OneOrMore)
    operon = RelationshipTo(Operon, HAS, cardinality=OneOrMore)
    gene = RelationshipTo(Gene, HAS, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, HAS, model=OperonRelationship)


class RegulatoryInteraction(Node):
    # properties
    regulator = StringProperty(required=True)
    operon = StringProperty(required=True)
    genes = ArrayProperty(StringProperty(), required=True)
    tfbss = ArrayProperty(StringProperty())
    effectors = ArrayProperty(StringProperty())
    regulatory_effect = StringProperty(required=True)

    # relationships
    source = RelationshipTo(Source, HAS, model=SourceRelationship)
    evidence = RelationshipTo(Evidence, HAS)
    publication = RelationshipTo(Publication, HAS)
    organism = RelationshipTo(Organism, HAS, cardinality=One)
    regulator_in = RelationshipTo(Regulator, HAS, cardinality=One)
    operon_out = RelationshipTo(Operon, HAS, cardinality=One)
    gene = RelationshipTo(Gene, HAS, cardinality=OneOrMore, model=OperonRelationship)
    tfbs = RelationshipTo(TFBS, HAS, model=OperonRelationship)
    effector = RelationshipTo(Effector, HAS)
