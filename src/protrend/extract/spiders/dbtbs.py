from typing import Union
from urllib import parse

from itemloaders import ItemLoader
from scrapy import Request, Spider
from scrapy.http import Response

from protrend.extract.items.dbtbs import TranscriptionFactorItem, GeneItem, TFBSItem
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

        tf_item = tf_loader.load_item()
        tf_loader = ItemLoader(item=tf_item)

        # regulatory interactions
        ri_table = response.xpath("//html/body/table[2]/tr")
        ri_header = ri_table[0]
        ri_body = ri_table[1:]

        cols = ri_header.xpath(".//th/a/text()").getall() + ri_header.xpath(".//th/text()").getall()

        if len(cols) > 6:
            mechanism = 'tf'
            operon_xpath = ".//td[1]/a/@href"
            gene_xpath = ".//td[2]/text()"
            tfbs_xpath = ".//td[7]/text()"
            absolute_position_xpath = ".//td[5]/text()"
            sequence_xpath = ".//td[7]/tt"
            pubmed_xpath = ".//td[8]/a/@href"

        else:
            mechanism = 'sigma'
            operon_xpath = ".//td[1]/a/@href"
            gene_xpath = ".//td[2]/text()"
            tfbs_xpath = ".//td[5]/text()"
            absolute_position_xpath = ".//td[3]/text()"
            sequence_xpath = ".//td[5]/tt"
            pubmed_xpath = ".//td[6]/a/@href"

        genes_items = []
        tfbs_items = []
        for ri in ri_body:

            gene = ri.xpath(gene_xpath).get()
            tfbs = ri.xpath(tfbs_xpath).get()

            if gene and tfbs:

                # gene
                gene_loader = ItemLoader(item=GeneItem(), selector=ri)
                gene_loader.add_xpath("name", gene_xpath)

                relative_url = ri.xpath(operon_xpath).get()
                url = parse.urljoin(response.url, relative_url)
                gene_loader.add_value("url", url)

                if mechanism == 'tf':
                    gene_loader.add_xpath("regulation", ".//td[4]/text()")
                else:
                    gene_loader.add_value("regulation", 'Promoter')

                gene_loader.add_xpath("pubmed", pubmed_xpath)

                # tfbs
                tfbs_loader = ItemLoader(item=TFBSItem(), selector=ri)

                tfbs_loader.add_value("identifier", tf_item.get("name"))
                tfbs_loader.add_xpath("identifier", gene_xpath)
                if mechanism == 'tf':
                    tfbs_loader.add_xpath("identifier", ".//td[4]/text()")
                else:
                    tfbs_loader.add_value("identifier", 'Promoter')

                tfbs_loader.add_xpath("identifier", absolute_position_xpath)

                relative_url = ri.xpath(operon_xpath).get()
                url = parse.urljoin(response.url, relative_url)
                tfbs_loader.add_value("url", url)

                if mechanism == 'tf':
                    tfbs_loader.add_xpath("regulation", ".//td[4]/text()")
                else:
                    tfbs_loader.add_value("regulation", 'Promoter')
                tfbs_loader.add_xpath("absolute_position", absolute_position_xpath)

                sequence_selector = ri.xpath(sequence_xpath)
                if sequence_selector:
                    sequence = sequence_selector.xpath("string()").extract()
                    if sequence:
                        tfbs_loader.add_value('sequence', sequence)

                tfbs_loader.add_xpath("pubmed", pubmed_xpath)

                # interaction
                gene_item = gene_loader.load_item()
                tfbs_item = tfbs_loader.load_item()

                gene_loader = ItemLoader(item=gene_item)
                tfbs_loader = ItemLoader(item=tfbs_item)

                gene_loader.add_value("tf", tf_item.get('name'))
                gene_loader.add_value("tfbs", tfbs_item.get('identifier'))

                tfbs_loader.add_value("tf", tf_item.get('name'))
                tfbs_loader.add_value("gene", gene_item.get('name'))

                gene_item = gene_loader.load_item()
                genes_items.append(gene_item)

                tfbs_item = tfbs_loader.load_item()
                tfbs_items.append(tfbs_item)

                tf_loader.add_value('gene', gene_item.get('name'))
                tfbs_loader.add_value('tfbs', tfbs_item.get('identifier'))

        tf_item = tf_loader.load_item()

        yield from genes_items
        yield from tfbs_items
        yield tf_item
