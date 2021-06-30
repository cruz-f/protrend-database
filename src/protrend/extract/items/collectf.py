from itemloaders.processors import TakeFirst, MapCompose, Join
from scrapy import Item, Field
from w3lib.html import remove_tags

from protrend.extract.processors import CollecTFProcessors


class TaxonomyItem(Item):
    taxonomy_id = Field(input_processor=MapCompose(CollecTFProcessors.process_tax_onclick),
                        output_processor=TakeFirst())
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
    operon_id = Field(output_processor=Join(separator='_'))

    # relationships
    regulon = Field()
    gene = Field()
    tfbs = Field()


class GeneItem(Item):
    locus_tag = Field(input_processor=MapCompose(str.strip),
                      output_processor=TakeFirst())

    # relationships
    regulon = Field()
    operon = Field()
    tfbs = Field()


class TFBSItem(Item):
    tfbs_id = Field(input_processor=MapCompose(CollecTFProcessors.process_site_identifier),
                    output_processor=Join(separator='_'))
    site_start = Field(output_processor=TakeFirst())
    site_end = Field(output_processor=TakeFirst())
    site_strand = Field(output_processor=TakeFirst())
    sequence = Field(output_processor=TakeFirst())
    mode = Field(input_processor=MapCompose(CollecTFProcessors.process_mode), output_processor=TakeFirst())
    pubmed = Field(input_processor=MapCompose(CollecTFProcessors.process_pubmed))

    # relationships
    organism = Field(output_processor=TakeFirst())
    regulon = Field()
    operon = Field()
    gene = Field()
    experimental_evidence = Field()


class ExperimentalEvidenceItem(Item):
    exp_id = Field(input_processor=MapCompose(str.strip, CollecTFProcessors.process_evidence_identifier),
                   output_processor=TakeFirst())

    # relationships
    regulon = Field()
    tfbs = Field()


class TranscriptionFactorItem(Item):
    name = Field(input_processor=MapCompose(str.strip), output_processor=TakeFirst())
    family = Field(input_processor=MapCompose(str.strip), output_processor=TakeFirst())
    description = Field(input_processor=MapCompose(remove_tags))
    pubmed = Field(input_processor=MapCompose(remove_tags, CollecTFProcessors.process_evidence_identifier))

    # relationships
    regulon = Field()
