import sys

from scrapy import cmdline

sys.path.insert(0, r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src')


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