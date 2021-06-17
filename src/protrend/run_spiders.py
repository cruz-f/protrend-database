from scrapy import cmdline


def run_spider(spider: str, log_file: str, urls: str = None):
    if urls:
        cmdline.execute(f"scrapy crawl {spider} -a urls={urls} --logfile {log_file}".split())

    else:
        cmdline.execute(f"scrapy crawl {spider}".split())


def run_regprecise(urls: str = None):
    log_file = 'regprecise.log'

    run_spider("regprecise", log_file, urls)


if __name__ == "__main__":
    regprecise_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                       "https://regprecise.lbl.gov/collections_tf.jsp",
                       "https://regprecise.lbl.gov/collections_tffam.jsp",
                       "https://regprecise.lbl.gov/collections_rfam.jsp",
                       "https://regprecise.lbl.gov/collections_effector.jsp",
                       "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_regprecise(regprecise_urls[0])
