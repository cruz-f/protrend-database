from itemloaders.processors import TakeFirst
from scrapy import Item, Field


class TaxonomyItem(Item):

    identifier = Field(output_processor=TakeFirst())
    name = Field(output_processor=TakeFirst())
    url = Field(output_processor=TakeFirst())

    # relationships
    organism = Field()


class OrganismItem(Item):

    name = Field(output_processor=TakeFirst())
    genome_accession = Field(output_processor=TakeFirst())

    # relationships
    taxonomy = Field(output_processor=TakeFirst())
    regulon = Field()
    tfbs = Field()


class RegulonItem(Item):

    uniprot_accession = Field(output_processor=TakeFirst())
    name = Field(output_processor=TakeFirst())
    url = Field(output_processor=TakeFirst())

    # relationships
    organism = Field(output_processor=TakeFirst())
    transcription_factor = Field(output_processor=TakeFirst())
    operon = Field()
    gene = Field()
    tfbs = Field()
    experimental_evidence = Field()


class OperonItem(Item):

    identifier = Field(output_processor=TakeFirst())
    name = Field(output_processor=TakeFirst())

    # relationships
    regulon = Field()
    gene = Field()
    tfbs = Field()


class GeneItem(Item):

    locus_tag = Field(output_processor=TakeFirst())

    # relationships
    regulon = Field()
    operon = Field()
    tfbs = Field()


class TFBSItem(Item):

    identifier = Field(output_processor=TakeFirst())
    site_start = Field(output_processor=TakeFirst())
    site_end = Field(output_processor=TakeFirst())
    site_strand = Field(output_processor=TakeFirst())
    sequence = Field(output_processor=TakeFirst())
    mode = Field(output_processor=TakeFirst())
    pubmed = Field()

    # relationships
    organism = Field(output_processor=TakeFirst())
    regulon = Field()
    operon = Field()
    gene = Field()
    experimental_evidence = Field()


class ExperimentalEvidenceItem(Item):

    identifier = Field(output_processor=TakeFirst())
    description = Field(output_processor=TakeFirst())

    # relationships
    regulon = Field()
    tfbs = Field()


class TranscriptionFactorFamilyItem(Item):

    identifier = Field(output_processor=TakeFirst())
    name = Field(output_processor=TakeFirst())
    description = Field(output_processor=TakeFirst())
    pubmed = Field()
    url = Field(output_processor=TakeFirst())

    # relationships
    transcription_factor = Field()


class TranscriptionFactorItem(Item):

    name = Field(output_processor=TakeFirst())
    description = Field(output_processor=TakeFirst())
    pubmed = Field()

    # relationships
    transcription_factor_family = Field(output_processor=TakeFirst())
    regulon = Field()