from typing import Dict

from neomodel import (StringProperty,
                      IntegerProperty,
                      RelationshipTo,
                      ArrayProperty)

from protrend.models.node import Node
from protrend.models.version import Version


# noinspection PyAbstractClass
class CollecTFVersion(Version):

    database = RelationshipTo('Database', 'VERSIONING')

    taxonomy = RelationshipTo('Taxonomy', 'VERSIONING')
    organism = RelationshipTo('Organism', 'VERSIONING')
    regulon = RelationshipTo('Regulon', 'VERSIONING')
    operon = RelationshipTo('Operon', 'VERSIONING')
    gene = RelationshipTo('Gene', 'VERSIONING')
    tfbs = RelationshipTo('TFBS', 'VERSIONING')
    experimentalevidence = RelationshipTo('ExperimentalEvidence', 'VERSIONING')

    transcriptionfactor = RelationshipTo('TranscriptionFactor', 'VERSIONING')
    transcriptionfactorfamily = RelationshipTo('TranscriptionFactorFamily', 'VERSIONING')

    @property
    def versioned_nodes(self) -> Dict[str, RelationshipTo]:

        return {'Database': self.database,
                'Taxonomy': self.taxonomy,
                'Organism': self.organism,
                'Regulon': self.regulon,
                'Operon': self.operon,
                'Gene': self.gene,
                'TFBS': self.tfbs,
                'ExperimentalEvidence': self.experimentalevidence,
                'TranscriptionFactor': self.transcriptionfactor,
                'TranscriptionFactorFamily': self.transcriptionfactorfamily}


class CollecTFNode(Node):

    __abstract_node__ = True

    version = RelationshipTo(CollecTFVersion, 'VERSIONING')


# noinspection PyAbstractClass
class Database(CollecTFNode):
    name = StringProperty(default='collectf')
    url = StringProperty(default='http://www.collectf.org/browse/browse/')
    doi = StringProperty(default='10.1093/nar/gkt1123')
    authors = ArrayProperty(StringProperty(), default=['Sefa Kılıç', 'Ivan Erill'])
    description = StringProperty(default='CollecTF: a database of experimentally validated transcription '
                                         'factor-binding sites in Bacteria')

    taxonomy = RelationshipTo('Taxonomy', 'HAS')
    transcription_factor_family = RelationshipTo('TranscriptionFactorFamily', 'HAS')


# noinspection PyAbstractClass
class Taxonomy(CollecTFNode):
    identifier = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    database = RelationshipTo(Database, 'IS_FROM')
    organism = RelationshipTo('Organism', 'HAS')


# noinspection PyAbstractClass
class Organism(CollecTFNode):
    name = StringProperty(required=True)
    genome_accession = StringProperty()

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    regulon = RelationshipTo('Regulon', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Regulon(CollecTFNode):
    uniprot_accession = StringProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    organism = RelationshipTo(Organism, 'IS_FROM')
    transcription_factor = RelationshipTo('TranscriptionFactor', 'IS_FROM')
    operon = RelationshipTo('Operon', 'HAS')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')
    experimental_evidence = RelationshipTo('ExperimentalEvidence', 'IS_FROM')


# noinspection PyAbstractClass
class Operon(CollecTFNode):
    identifier = StringProperty(required=True)
    name = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Gene(CollecTFNode):
    locus_tag = StringProperty(required=True)

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class TFBS(CollecTFNode):
    identifier = StringProperty(required=True)
    site_start = IntegerProperty(required=True)
    site_end = IntegerProperty(required=True)
    site_strand = IntegerProperty(required=True)
    sequence = StringProperty(required=True)
    mode = StringProperty(required=True)
    pubmed = ArrayProperty(StringProperty())

    organism = RelationshipTo(Organism, 'IS_FROM')
    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    gene = RelationshipTo(Gene, 'IS_FROM')
    experimental_evidence = RelationshipTo('ExperimentalEvidence', 'IS_FROM')


# noinspection PyAbstractClass
class ExperimentalEvidence(CollecTFNode):
    identifier = StringProperty(required=True)
    description = StringProperty()


# noinspection PyAbstractClass
class TranscriptionFactorFamily(CollecTFNode):
    identifier = IntegerProperty(required=True)
    name = StringProperty(required=True)
    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())
    url = StringProperty()

    database = RelationshipTo(Database, 'IS_FROM')
    transcription_factor = RelationshipTo('TranscriptionFactor', 'HAS')


# noinspection PyAbstractClass
class TranscriptionFactor(CollecTFNode):
    name = StringProperty(required=True)
    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    transcription_factor_family = RelationshipTo('TranscriptionFactorFamily', 'IS_FROM')
    regulon = RelationshipTo(Regulon, 'HAS')
