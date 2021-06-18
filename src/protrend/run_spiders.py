from scrapy import cmdline


def run_spider(spider: str,
               log_file: str,
               urls: str = None,
               clear_version: bool = False,
               clear_sa: bool = False,
               clear_sa_schema: bool = False):

    if urls:
        cmdline.execute(f"scrapy crawl {spider} "
                        f"-a urls={urls} "
                        f"--logfile {log_file} "
                        f"-s clear_version={clear_version} "
                        f"-s clear_sa={clear_sa} "
                        f"-s clear_sa_schema={clear_sa_schema} ".split())

    else:
        cmdline.execute(f"scrapy crawl {spider} "
                        f"--logfile {log_file} "
                        f"-s clear_version={clear_version} "
                        f"-s clear_sa={clear_sa} "
                        f"-s clear_sa_schema={clear_sa_schema} ".split())


def run_reg_precise(urls: str = None,
                    clear_version: bool = False,
                    clear_sa: bool = False,
                    clear_sa_schema: bool = False):

    log_file = 'regprecise.log'

    run_spider("regprecise",
               log_file,
               urls,
               clear_version,
               clear_sa,
               clear_sa_schema)


if __name__ == "__main__":
    reg_precise_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                        "https://regprecise.lbl.gov/collections_tf.jsp",
                        "https://regprecise.lbl.gov/collections_tffam.jsp",
                        "https://regprecise.lbl.gov/collections_rfam.jsp",
                        "https://regprecise.lbl.gov/collections_effector.jsp",
                        "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_reg_precise(reg_precise_urls[0],
                    clear_sa_schema=True)
