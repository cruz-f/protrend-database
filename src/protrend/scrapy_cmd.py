from scrapy import cmdline


def run_spider(spider: str):

    cmdline.execute(f"scrapy crawl {spider}".split())


if __name__ == "__main__":

    pass