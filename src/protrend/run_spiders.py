from scrapy import cmdline


def run_spider(spider: str, urls: str = None):
    if urls:
        cmdline.execute(f"scrapy crawl {spider} -a urls={urls}".split())

    else:
        cmdline.execute(f"scrapy crawl {spider}".split())


def run_regprecise(urls: str = None):
    run_spider("regprecise", urls)


if __name__ == "__main__":
    urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
            "https://regprecise.lbl.gov/collections_tf.jsp",
            "https://regprecise.lbl.gov/collections_tffam.jsp",
            "https://regprecise.lbl.gov/collections_rfam.jsp",
            "https://regprecise.lbl.gov/collections_effector.jsp",
            "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_regprecise(urls[5])
