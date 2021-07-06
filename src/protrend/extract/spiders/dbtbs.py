from typing import Union

from itemloaders import ItemLoader
from scrapy import Request, Spider
from scrapy.http import Response

from protrend.extract.items.dbtbs import TranscriptionFactorItem, OperonItem, GeneItem, TFBSItem
from protrend.extract.utils import parsing_spider_arguments


class DBTBSSpider(Spider):
    name = "dbtbs"

    custom_settings = {
        'ITEM_PIPELINES': {'extract.pipelines.dbtbs.DBTBSPipeline': 800}
    }

    start_urls = ("https://dbtbs.hgc.jp/tfactable.html", )

    allowed_domains = ["dbtbs.hgc.jp"]

    def __init__(self, urls: Union[str, tuple, list, set] = None, *args, **kwargs):

        super(DBTBSSpider, self).__init__(*args, **kwargs)

        urls = parsing_spider_arguments(urls)

        if urls:
            self.start_urls = urls

    def start_requests(self):

        for url in self.start_urls:

            if url.endswith("tfactable.html"):
                yield Request(url, callback=self.parse_tfs)

    def parse_tfs(self, response: Response, **kwargs):

        domains_xpath = "//html/body/table[1]/tr"
        family_xpath = ".//th/a/text()"
        domain_xpath = ".//td[1]/a/text()"
        domain_description_xpath = ".//td[2]/text()"

        family = None
        domains = []

        domains_selector = response.xpath(domains_xpath)

        for domain_selector in domains_selector:

            current_family = domain_selector.xpath(family_xpath).get(None)

            if current_family:
                family = current_family

            else:
                domain = domain_selector.xpath(domain_xpath).get()
                domain_description = domain_selector.xpath(domain_description_xpath).get()
                domains.append((domain, domain_description, family))

        domains_tables_xpath = "//html/body/table"
        domains_tables = response.xpath(domains_tables_xpath)

        for (domain, domain_description, family), domain_table in zip(domains, domains_tables[1:]):

            tfs_xpath = ".//tr"
            tfs = domain_table.xpath(tfs_xpath)

            for tf in tfs:

                name_xpath = ".//td[1]/a/text()"
                name = tf.xpath(name_xpath)

                if name:
                    tf_loader = ItemLoader(item=TranscriptionFactorItem(), selector=tf)

                    tf_loader.add_xpath("name", name_xpath)

                    tf_loader.add_value("family", family)
                    tf_loader.add_value("domain", domain)
                    tf_loader.add_value("domain_description", domain_description)

                    description_xpath = ".//td[4]/text()"
                    tf_loader.add_xpath("description", description_xpath)

                    tf_item = tf_loader.load_item()

                    tf_link = tf.xpath(".//td[1]/a")[0]
                    cb_kwargs = dict(tf_item=tf_item)

                    yield response.follow(tf_link,
                                          callback=self.parse_tf,
                                          cb_kwargs=cb_kwargs)

    def parse_tf(self, response: Response, tf_item: TranscriptionFactorItem):

        tf_loader = ItemLoader(item=tf_item, selector=response)

        tf_loader.add_value("url", response.url)

        type_xpath = "//html/body/table[1]/tr[1]/td[1]/text()"
        tf_loader.add_xpath("type", type_xpath)

        swiss_xpath = "//html/body/table[1]/tr[2]/td[1]/text()"
        tf_loader.add_xpath("swiss_prot", swiss_xpath)

        subti_xpath = "//html/body/table[1]/tr[3]/td[1]/a/text()"
        tf_loader.add_xpath("subti_list", subti_xpath)

        consensus_xpath = "//html/body/table[1]/tr[4]/td[1]/text()"
        tf_loader.add_xpath("consensus_sequence", consensus_xpath)

        comment_xpath = "//html/body/table[1]/tr[5]/td[1]/text()"
        tf_loader.add_xpath("comment", comment_xpath)

        tf_loader.load_item()

        operons_xpath = "//html/body/table[2]/tr/td[1]/a"
        operons = response.xpath(operons_xpath)

        for operon in operons:
            operon_loader = ItemLoader(item=OperonItem(), selector=operon)

            name_xpath = ".//text()"
            operon_loader.add_xpath("name", name_xpath)

            operon_loader.add_value("tf", tf_item.get("name"))
            operon_item = operon_loader.load_item()

            tf_loader = ItemLoader(item=tf_item)
            tf_loader.add_value("operon", operon_item.get("name"))
            tf_loader.load_item()

            cb_kwargs = dict(tf_name=tf_item.get("name"), operon_item=operon_item)

            yield response.follow(operon,
                                  callback=self.parse_operon,
                                  cb_kwargs=cb_kwargs)

        yield tf_item

    def parse_operon(self, response: Response, tf_name: str, operon_item: OperonItem):

        operon_loader = ItemLoader(item=operon_item, selector=response)
        operon_loader.add_value("url", response.url)

        evidence_xpath = "//html/body/table[3]/tr[1]/td/text()"
        operon_loader.add_xpath("evidence", evidence_xpath)

        pubmed_xpath = "//html/body/table[3]/tr[2]/td/a/@href"
        operon_loader.add_xpath("pubmed", pubmed_xpath)

        comment_xpath = "//html/body/table[3]/tr[3]/td/text()"
        operon_loader.add_xpath("comment", comment_xpath)

        genes_items = self.parse_genes(response)
        tfbs_items = self.parse_tfbs(response, operon_item)

        for gene_item in genes_items:

            operon_loader.add_value("gene", gene_item.get("name"))

            gene_loader = ItemLoader(item=gene_item)
            gene_loader.add_value("tf", tf_name)
            gene_loader.add_value("operon", operon_item.get("name"))
            gene_loader.add_value("tfbs", [tfbs_item.get('identifier') for tfbs_item in tfbs_items])
            gene_loader.load_item()

        for tfbs_item in tfbs_items:
            operon_loader.add_value("tfbs", tfbs_item.get("identifier"))

            tfbs_loader = ItemLoader(item=tfbs_item)
            tfbs_loader.add_value("tf", tf_name)
            tfbs_loader.add_value("operon", operon_item.get("name"))
            tfbs_loader.add_value("gene", [gene_item.get('name') for gene_item in genes_items])
            tfbs_loader.load_item()

        yield operon_loader.load_item()
        yield from genes_items
        yield from tfbs_items

    @staticmethod
    def parse_genes(response: Response):

        genes_xpath = "//html/body/table[2]/tr"
        genes = response.xpath(genes_xpath)

        genes_items = []

        for gene in genes:

            name_xpath = ".//td[1]/text()"
            name = gene.xpath(name_xpath)

            if name:

                gene_loader = ItemLoader(item=GeneItem(), selector=gene)

                gene_loader.add_xpath("name", name_xpath)
                gene_loader.add_value("url", response.url)

                synonyms_xpath = ".//td[2]/text()"
                gene_loader.add_xpath("synonyms", synonyms_xpath)

                strand_xpath = ".//td[3]/text()"
                gene_loader.add_xpath("strand", strand_xpath)

                position_xpath = ".//td[4]/text()"
                gene_loader.add_xpath("position", position_xpath)

                function_xpath = ".//td[5]/text()"
                gene_loader.add_xpath("function", function_xpath)

                cog_id_xpath = ".//td[6]/a/text()"
                gene_loader.add_xpath("cog_id", cog_id_xpath)

                conversed_groups_xpath = ".//td[7]/text()"
                gene_loader.add_xpath("conversed_groups", conversed_groups_xpath)

                gene_item = gene_loader.load_item()
                genes_items.append(gene_item)

        return genes_items

    @staticmethod
    def parse_tfbs(response: Response, operon_item: OperonItem):

        tfbs_xpath = "//html/body/table[4]/tr"
        tfbs = response.xpath(tfbs_xpath)

        tfbs_items = []

        for bs in tfbs:

            tf_xpath = ".//td[1]/a/text()"
            tf = bs.xpath(tf_xpath)

            if tf:
                tfbs_loader = ItemLoader(item=TFBSItem(), selector=bs)

                tfbs_loader.add_xpath("identifier", tf_xpath)
                tfbs_loader.add_value("url", response.url)
                tfbs_loader.add_value("identifier", operon_item.get("name"))

                regulation_xpath = ".//td[2]/text()"
                tfbs_loader.add_xpath("identifier", regulation_xpath)
                tfbs_loader.add_xpath("regulation", regulation_xpath)

                location_xpath = ".//td[3]/tt/text()"
                tfbs_loader.add_xpath("location", location_xpath)

                absolute_position_xpath = ".//td[4]/text()"
                tfbs_loader.add_xpath("identifier", absolute_position_xpath)
                tfbs_loader.add_xpath("absolute_position", absolute_position_xpath)

                sequence_xpath = ".//td[5]/tt/text()"
                tfbs_loader.add_xpath("sequence", sequence_xpath)

                pubmed_xpath = ".//td[6]/a/@href"
                tfbs_loader.add_xpath("pubmed", pubmed_xpath)

                tfbs_item = tfbs_loader.load_item()
                tfbs_items.append(tfbs_item)

        return tfbs_items
