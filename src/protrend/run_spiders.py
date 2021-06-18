from scrapy import cmdline


def run_spider(spider: str,
               log_file: str,
               urls: str = None,
               version: str = None,
               clear_version: bool = False,
               clear_db: bool = False,
               clear_schema: bool = False):

    if urls:

        if version is None:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"-a urls={urls} "
                            f"--logfile {log_file} "
                            f"-s clear_version={clear_version} "
                            f"-s clear_db={clear_db} "
                            f"-s clear_schema={clear_schema} ".split())

        else:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"-a urls={urls} "
                            f"--logfile {log_file} "
                            f"-s version={version} "
                            f"-s clear_version={clear_version} "
                            f"-s clear_db={clear_db} "
                            f"-s clear_schema={clear_schema} ".split())

    else:
        if version is None:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"--logfile {log_file} "
                            f"-s clear_version={clear_version} "
                            f"-s clear_db={clear_db} "
                            f"-s clear_schema={clear_schema} ".split())

        else:
            cmdline.execute(f"scrapy crawl {spider} "
                            f"--logfile {log_file} "
                            f"-s version={version} "
                            f"-s clear_version={clear_version} "
                            f"-s clear_db={clear_db} "
                            f"-s clear_schema={clear_schema} ".split())


def run_reg_precise(urls: str = None,
                    version: str = None,
                    clear_version: bool = False,
                    clear_db: bool = False,
                    clear_schema: bool = False):

    log_file = 'regprecise.log'

    run_spider(spider="regprecise",
               log_file=log_file,
               urls=urls,
               version=version,
               clear_version=clear_version,
               clear_db=clear_db,
               clear_schema=clear_schema)


if __name__ == "__main__":
    reg_precise_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                        "https://regprecise.lbl.gov/collections_tf.jsp",
                        "https://regprecise.lbl.gov/collections_tffam.jsp",
                        "https://regprecise.lbl.gov/collections_rfam.jsp",
                        "https://regprecise.lbl.gov/collections_effector.jsp",
                        "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_reg_precise(version='0.0.0', clear_schema=True)
