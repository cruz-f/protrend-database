from collections import namedtuple

from scrapy import Spider, Request
from scrapy.http import Response
from scrapy.loader import ItemLoader

from protrend.extract.utils import parsing_spider_arguments
from protrend.extract.items import TaxonomyItem, GenomeItem, RegulonItem, OperonItem, TFBSItem, GeneItem, \
    TranscriptionFactorItem, RegulogItem

OperonGeneTFBS = namedtuple("OperonGeneTFBS",
                            ["operon", "genes", "tfbs"],
                            defaults=[None, None, None])


class RegPreciseSpider(Spider):
    name = "regprecise"

    start_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                  "https://regprecise.lbl.gov/collections_tf.jsp",
                  "https://regprecise.lbl.gov/collections_tffam.jsp",
                  "https://regprecise.lbl.gov/collections_rfam.jsp",
                  "https://regprecise.lbl.gov/collections_effector.jsp",
                  "https://regprecise.lbl.gov/collections_pathway.jsp")

    allowed_domains = ["regprecise.lbl.gov"]

    def __init__(self, urls=None, *args, **kwargs):

        super(RegPreciseSpider, self).__init__(*args, **kwargs)

        urls = parsing_spider_arguments(urls)

        if urls:
            self.start_urls = urls

    def start_requests(self):

        for url in self.start_urls:

            if url.endswith("collections_tax.jsp"):
                yield Request(url, self.parse_collections_tax)

            elif url.endswith("collections_tf.jsp"):
                yield Request(url, self.parse_collections_tf)

            elif url.endswith("collections_tffam.jsp"):
                yield Request(url, self.parse_collections_tffam)

            elif url.endswith("collections_rfam.jsp"):
                yield Request(url, self.parse_collections_rfam)

            elif url.endswith("collections_effector.jsp"):
                yield Request(url, self.parse_collections_effector)

            elif url.endswith("collections_pathway.jsp"):
                yield Request(url, self.parse_collections_pathway)

            else:
                yield Request(url, self.parse)

    def parse_collections_tax(self, response: Response):

        taxa_xpath = "//*[@id='content']/div/div[3]/table/tbody/tr/td[1]/a"
        taxa = response.xpath(taxa_xpath)

        for taxon in taxa:
            taxon_loader = ItemLoader(item=TaxonomyItem(),
                                      selector=taxon)

            anchor_href_xpath = ".//@href"
            taxon_loader.add_xpath("collection_id", anchor_href_xpath)

            anchor_text_xpath = ".//text()"
            taxon_loader.add_xpath("name", anchor_text_xpath)

            taxon_item = taxon_loader.load_item()
            cb_kwargs = dict(taxon_item=taxon_item)
            yield response.follow(taxon,
                                  callback=self.parse_genome,
                                  cb_kwargs=cb_kwargs)

            yield taxon_item

    def parse_genome(self, response: Response, taxon_item: TaxonomyItem):

        taxon_loader = ItemLoader(item=taxon_item,
                                  response=response)

        taxon_loader.add_value("url", response.url)

        genomes_xpath = "//*[@id='content']/div/table/tbody/tr/td[1]/a"
        genomes = response.xpath(genomes_xpath)

        for genome in genomes:
            genome_loader = ItemLoader(item=GenomeItem(),
                                       selector=genome)

            genome_loader.add_value("taxonomy", taxon_item.collection_id)

            anchor_href_xpath = ".//@href"
            genome_loader.add_xpath("genome_id", anchor_href_xpath)
            taxon_loader.add_xpath("genome", anchor_href_xpath)

            anchor_text_xpath = ".//text()"
            genome_loader.add_xpath("name", anchor_text_xpath)

            genome_item = genome_loader.load_item()
            cb_kwargs = dict(genome_item=genome_item)

            yield response.follow(genome,
                                  callback=self.parse_regulon,
                                  cb_kwargs=cb_kwargs)

            yield genome_item

    def parse_regulon(self, response: Response, genome_item: GenomeItem):

        genome_loader = ItemLoader(item=genome_item,
                                   response=response)

        genome_loader.add_value("url", response.url)

        regulons_xpath = "//*[@id='content']/div[2]/div/table/tbody/tr/td[2]/a"
        regulons = response.xpath(regulons_xpath)

        for regulon in regulons:
            regulon_loader = ItemLoader(item=RegulonItem(),
                                        selector=regulon)

            regulon_loader.add_value("genome", genome_item.genome_id)

            anchor_href_xpath = ".//@href"
            regulon_loader.add_xpath("regulon_id", anchor_href_xpath)
            genome_loader.add_xpath("regulon", anchor_href_xpath)

            anchor_text_xpath = ".//text()"
            regulon_loader.add_xpath("name", anchor_text_xpath)

            regulon_item = regulon_loader.load_item()
            cb_kwargs = dict(regulon_item=regulon_item)

            yield response.follow(regulon,
                                  callback=self.parse_regulon_page,
                                  cb_kwargs=cb_kwargs)

            yield regulon_item

    def parse_regulon_page(self, response: Response, regulon_item: RegulonItem):

        regulon_loader = ItemLoader(item=regulon_item,
                                    response=response)

        regulon_loader.add_value("url", response.url)

        regulon_item = regulon_loader.load_item()

        self.parse_regulon_properties(response=response,
                                      regulon_item=regulon_item)

        self.parse_regulon_collections(response=response, regulon_item=regulon_item)

        operon_genes_tfbs = self.parse_operons_genes_tfbs(response=response, regulon_item=regulon_item)

        for (operon, genes, tfbs) in operon_genes_tfbs:
            yield from tfbs

            yield from genes

            yield operon

    @staticmethod
    def parse_regulon_properties(response: Response, regulon_item: RegulonItem):

        regulon_loader = ItemLoader(item=regulon_item,
                                    response=response)

        regulator_type = "//*[@id='propblock']/table/tbody/tr[1]/td[2]//text()"
        regulon_loader.add_xpath("regulator_type", regulator_type)

        regulator_locus_tag = "//*[@id='propblock']/table/tbody/tr[2]/td[2]/a//text()"
        regulon_loader.add_xpath("regulator_locus_tag", regulator_locus_tag)

        regulator_family = "//*[@id='propblock']/table/tbody/tr[3]/td[2]//text()"
        regulon_loader.add_xpath("regulator_family", regulator_family)

        regulation_mode = ", '//*[@id='propblock']/table/tbody/tr[4]/td[2]//text()"
        regulon_loader.add_xpath("regulation_mode", regulation_mode)

        biological_process = "//*[@id='propblock']/table/tbody/tr[5]/td[2]//text()"
        regulon_loader.add_xpath("biological_process", biological_process)

        regulation_effector = "//*[@id='propblock']/table/tbody/tr[6]/td[2]//text()"
        regulon_loader.add_xpath("regulation_effector", regulation_effector)

        regulation_regulog = "//*[@id='propblock']/table/tbody/tr[7]/td[2]/a//text()"
        regulon_loader.add_xpath("regulation_regulog", regulation_regulog)

        regulog = "//*[@id='propblock']/table/tbody/tr[7]/td[2]/a//@href"
        regulon_loader.add_xpath("regulog", regulog)

        regulon_loader.load_item()

    @staticmethod
    def parse_regulon_collections(response: Response, regulon_item: RegulonItem):

        regulon_loader = ItemLoader(item=regulon_item,
                                    response=response)

        collections_xpath = "//*[@id='content']/div[6]/ul/li"
        collections = response.xpath(collections_xpath)

        for collection in collections:

            collection_text = collection.xpath(".//text()").get()

            if collection_text:

                if "taxonomy" in collection_text:

                    regulon_loader.add_value("taxonomy", ".//@href")

                elif "trascription factor" in collection_text:

                    regulon_loader.add_value("transcription_factor", ".//@href")

                elif "TF family" in collection_text:

                    regulon_loader.add_value("tf_family", ".//@href")

                elif "RNA motif" in collection_text:

                    regulon_loader.add_value("rna_family", ".//@href")

                elif "effector" in collection_text:

                    regulon_loader.add_value("effector", ".//@href")

                elif "pathway" in collection_text:

                    regulon_loader.add_value("pathway", ".//@href")

        regulon_loader.load_item()

    @staticmethod
    def parse_operons_genes_tfbs(response: Response, regulon_item: RegulonItem):

        operon_xpath = "//*[@id='operontbl']/tbody/tr/td/div[contains(@class, 'operon')]"
        operons = response.xpath(operon_xpath)

        operon_gene_tfbs = []

        # loading main attributes and creating operon identifier
        for operon in operons:

            operon_loader = ItemLoader(item=OperonItem())

            operon_tfbs = []
            operon_genes = []

            tooltips_xpath = ".//div[contains(@class, 'tooltip')]"
            tooltips = operon.xpath(tooltips_xpath)

            for tooltip in tooltips:

                tfbs_item = RegPreciseSpider.parse_tfbs(tooltip=tooltip)
                if tfbs_item:
                    operon_tfbs.append(tfbs_item)

                gene_item = RegPreciseSpider.parse_gene(tooltip=tooltip)

                if gene_item:
                    operon_genes.append(gene_item)

                    # operon identifier is created according to the operon genes
                    operon_loader.add_value("operon_id", gene_item.locus_tag)
                    operon_loader.add_value("name", gene_item.name)

            operon_item = operon_loader.load_item()
            operon_gene_tfbs.append(OperonGeneTFBS(operon=operon_item,
                                                   genes=operon_genes,
                                                   tfbs=operon_tfbs))

        # setting up connections and creating tfbs identifier
        regulon_loader = ItemLoader(item=regulon_item)

        for (operon, genes, tfbs) in operon_gene_tfbs:

            regulon_loader.add_value("operon", operon.operon_id)

            operon_loader = ItemLoader(item=operon)

            operon_loader.add_value("url", response.url)
            operon_loader.add_value("regulon", regulon_item.regulon_id)

            for tfbs_item in tfbs:

                tfbs_loader = ItemLoader(item=tfbs_item)

                tfbs_loader.add_value("tfbs_id", tfbs_item.position)
                tfbs_loader.add_value("tfbs_id", operon.operon_id)

                tfbs_loader.add_value("url", response.url)
                tfbs_loader.add_value("regulon", regulon_item.regulon_id)
                tfbs_loader.add_value("operon", operon.operon_id)

                for gene_item in genes:
                    tfbs_loader.add_value("gene", gene_item.locus_tag)

                tfbs_loader.load_item()

                operon_loader.add_value("tfbs", tfbs_item.tfbs_id)

                regulon_loader.add_value("tfbs", tfbs_item.tfbs_id)

            for gene_item in genes:

                gene_loader = ItemLoader(item=gene_item)

                gene_loader.add_value("url", response.url)
                gene_loader.add_value("regulon", regulon_item.regulon_id)
                gene_loader.add_value("operon", operon.operon_id)

                for tfbs_item in tfbs:
                    gene_loader.add_value("tfbs", tfbs_item.tfbs_id)

                operon_loader.add_value("gene", gene_item.locus_tag)

                regulon_loader.add_value("gene", gene_item.locus_tag)

        return operon_gene_tfbs

    @staticmethod
    def parse_tfbs(tooltip):

        tooltip_site_img = tooltip.xpath(".//div[contains(@class, 'site_img')]")

        if tooltip_site_img:

            tfbs_props = tooltip.xpath(".//span//text()")

            if tfbs_props:
                tfbs_loader = ItemLoader(item=TFBSItem())

                position, score, sequence = tfbs_props.getall()
                tfbs_loader.add_value("position", position)
                tfbs_loader.add_value("score", score)
                tfbs_loader.add_value("sequence", sequence)

                return tfbs_loader.load_item()

        return

    @staticmethod
    def parse_gene(tooltip):

        tooltip_gene_img = tooltip.xpath(".//div[contains(@class, 'gene_img')]")

        if tooltip_gene_img:

            genes_props = tooltip.xpath(".//span//text()")

            if genes_props:
                gene_loader = ItemLoader(item=GeneItem())

                locus_tag, name, func = genes_props.getall()
                gene_loader.add_value("locus_tag", locus_tag)
                gene_loader.add_value("name", name)
                gene_loader.add_value("function", func)

                return gene_loader.load_item()

        return

    # TODO: missing checking parsing of tfs
    def parse_collections_tf(self, response):

        # TODO: check valid xpaths
        tf_xpath = "//*[@id='content']/div/div[3]/table/tbody/tr/td[1]/a"
        tfs = response.xpath(tf_xpath)

        for tf in tfs:
            tf_loader = ItemLoader(item=TranscriptionFactorItem(),
                                   selector=tf)

            anchor_href_xpath = ".//@href"
            tf_loader.add_xpath("collection_id", anchor_href_xpath)

            anchor_text_xpath = ".//text()"
            tf_loader.add_xpath("name", anchor_text_xpath)

            tf_item = tf_loader.load_item()
            cb_kwargs = dict(tf_item=tf_item)
            yield response.follow(tf,
                                  callback=self.parse_regulog,
                                  cb_kwargs=cb_kwargs)

            yield tf_item

    def parse_regulog(self, response, tf_item):

        tf_loader = ItemLoader(item=tf_item,
                               response=response)

        tf_loader.add_value("url", response.url)

        # TODO: parsing description, pubmed

        # TODO: check valid xpaths
        regulogs_xpath = "//*[@id='content']/div/table/tbody/tr/td[1]/a"
        regulogs = response.xpath(regulogs_xpath)

        for regulog in regulogs:
            regulog_loader = ItemLoader(item=RegulogItem(),
                                        selector=regulog)

            regulog_loader.add_value("transcription_factor", tf_item.collection_id)

            anchor_href_xpath = ".//@href"
            regulog_loader.add_xpath("regulog_id", anchor_href_xpath)
            tf_loader.add_xpath("regulog", anchor_href_xpath)

            anchor_text_xpath = ".//text()"
            regulog_loader.add_xpath("name", anchor_text_xpath)

            regulog_item = regulog_loader.load_item()
            cb_kwargs = dict(regulog_item=regulog_item)

            yield response.follow(regulog,
                                  callback=self.parse_regulog_page,
                                  cb_kwargs=cb_kwargs)

            yield regulog_item

    def parse_regulog_page(self, response, regulog_item):

        regulog_loader = ItemLoader(item=regulog_item,
                                    response=response)

        regulog_loader.add_value("url", response.url)

        regulog_item = regulog_loader.load_item()

        self.parse_regulog_properties(response=response,
                                      regulog_item=regulog_item)

        self.parse_regulog_collections(response=response, regulog_item=regulog_item)

        self.parse_regulog_regulons(response=response, regulog_item=regulog_item)

    @staticmethod
    def parse_regulog_properties(response, regulog_item):

        # TODO: check valid xpaths

        regulog_loader = ItemLoader(item=regulog_item,
                                    response=response)

        regulator_type = "//*[@id='propblock']/table/tbody/tr[1]/td[2]//text()"
        regulog_loader.add_xpath("regulator_type", regulator_type)

        regulator_family = "//*[@id='propblock']/table/tbody/tr[2]/td[2]//text()"
        regulog_loader.add_xpath("regulator_family", regulator_family)

        regulation_mode = ", '//*[@id='propblock']/table/tbody/tr[3]/td[2]//text()"
        regulog_loader.add_xpath("regulation_mode", regulation_mode)

        biological_process = "//*[@id='propblock']/table/tbody/tr[4]/td[2]//text()"
        regulog_loader.add_xpath("biological_process", biological_process)

        regulation_effector = "//*[@id='propblock']/table/tbody/tr[5]/td[2]//text()"
        regulog_loader.add_xpath("regulation_effector", regulation_effector)

        phylum = "//*[@id='propblock']/table/tbody/tr[6]/td[2]/a//text()"
        regulog_loader.add_xpath("phylum", phylum)

        regulog_loader.load_item()

    def parse_regulog_collections(self, response, regulog_item):
        # TODO: check valid xpaths
        return self.parse_regulon_collections(response, regulog_item)

    def parse_regulog_regulons(self, response, regulog_item):

        # TODO: missing parsing of regulog regulons
        pass

    def parse_collections_tffam(self, response, **kwargs):
        # TODO: missing parsing of tffam
        pass

    def parse_collections_rfam(self, response, **kwargs):
        # TODO: missing parsing of rfam
        pass

    def parse_collections_effector(self, response, **kwargs):
        # TODO: missing parsing of effector
        pass

    def parse_collections_pathway(self, response, **kwargs):
        # TODO: missing parsing of pathway
        pass

    def parse(self, response, **kwargs):

        yield
