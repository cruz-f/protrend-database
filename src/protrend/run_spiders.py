from scrapy import cmdline


def run_spider(spider: str,
               logfile: str,
               user_name: str,
               password: str,
               ip: str,
               port: str,
               db_name: str,
               dbms: str,
               import_folder: str,
               version: str = None,
               urls: str = None):
    """


    Crawler settings
    :param spider:
    :param logfile:

    DB settings
    :param user_name:
    :param password:
    :param ip:
    :param port:
    :param db_name:
    :param dbms:

    DB pipeline
    :param import_folder:
    :param version:

    Spider settings
    :param urls:
    :return:
    """

    arguments = ["scrapy",
                 "crawl",
                 spider,
                 "-s", f"LOG_FILE={logfile}",
                 "-s", f"user_name={user_name}",
                 "-s", f"password={password}",
                 "-s", f"ip={ip}",
                 "-s", f"port={port}",
                 "-s", f"db_name={db_name}",
                 "-s", f"dbms={dbms}",
                 "-s", f"import_folder={import_folder}"]

    if version:
        arguments.append("-s")
        arguments.append(f"version={version}")

    if urls:
        arguments.append("-a")
        arguments.append(f"urls={urls}")

    return cmdline.execute(arguments)


if __name__ == "__main__":
    reg_precise_urls = ("https://regprecise.lbl.gov/collections_tax.jsp",
                        "https://regprecise.lbl.gov/collections_tf.jsp",
                        "https://regprecise.lbl.gov/collections_tffam.jsp",
                        "https://regprecise.lbl.gov/collections_rfam.jsp",
                        "https://regprecise.lbl.gov/collections_effector.jsp",
                        "https://regprecise.lbl.gov/collections_pathway.jsp")

    run_spider(spider='regprecise',
               logfile='regprecise.log',
               user_name='neo4j',
               password='regprecise',
               ip='localhost',
               port='7687',
               db_name='neo4j',
               dbms=r'C:\Users\BiSBII\.Neo4jDesktop\relate-data\dbmss\dbms-0d680ec0-bb15-4f6e-9992-aa7f72201baf',
               import_folder=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src\protrend\extract\import\regprecise',
               version='0.0.0',
               urls=reg_precise_urls[3])
