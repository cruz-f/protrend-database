import scrapy

from protrend.extract.extract_utils import parsing_spider_arguments


class RegPreciseSpider(scrapy.Spider):
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
            yield scrapy.Request(url, self.parse)

    def parse_collections_tax(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse_collections_tf(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse_collections_tffam(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse_collections_rfam(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse_collections_effector(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse_collections_pathway(self, response, **kwargs):
        # //*[(@id = "content")]//a

        pass

    def parse(self, response, **kwargs):
        for quote in response.css('div.quote'):
            yield {
                'text': quote.css('span.text::text').get(),
                'author': quote.css('small.author::text').get(),
            }

        next_page = response.css('li.next a::attr(href)').get()
        if next_page is not None:
            yield response.follow(next_page, self.parse)
