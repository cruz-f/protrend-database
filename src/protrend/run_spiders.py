from scrapy import cmdline


def run_spider(spider: str,
               log_file: str,
               urls: str = None,
               version: str = None):

    if urls:

        if version is None:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"-a urls={urls} "
                            f"--logfile {log_file} ".split())

        else:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"-a urls={urls} "
                            f"--logfile {log_file} "
                            f"-s version={version} ".split())

    else:
        if version is None:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"--logfile {log_file} ".split())

        else:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"--logfile {log_file} "
                            f"-s version={version} ".split())


def run_reg_precise(urls: str = None,
                    version: str = None):

    log_file = 'regprecise.log'

    run_spider(spider="regprecise",
               log_file=log_file,
               urls=urls,
               version=version)


if __name__ == "__main__":
    reg_precise_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                        "https://regprecise.lbl.gov/collections_tf.jsp",
                        "https://regprecise.lbl.gov/collections_tffam.jsp",
                        "https://regprecise.lbl.gov/collections_rfam.jsp",
                        "https://regprecise.lbl.gov/collections_effector.jsp",
                        "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_reg_precise(version='0.0.0')
