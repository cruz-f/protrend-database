import os
import shutil

from scrapy import cmdline


def run_spider(spider: str,
               staging_area: str,
               version: str = None,
               urls: str = None):
    """


    Crawler settings
    :param spider:

    Spider pipeline
    :param staging_area:
    :param version:

    Spider settings
    :param urls:
    :return:
    """

    if not os.path.exists(staging_area):
        os.makedirs(staging_area)

    spider_sa = os.path.join(staging_area, spider)

    if not os.path.exists(spider_sa):
        os.makedirs(spider_sa)

    spider_sa_version = os.path.join(staging_area, spider, version)

    if not os.path.exists(spider_sa_version):
        os.makedirs(spider_sa_version)

    else:
        for filename in os.listdir(spider_sa_version):
            file_path = os.path.join(spider_sa_version, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')

    logfile = fr'{spider_sa_version}\{spider}.log'

    try:
        if os.path.isfile(logfile) or os.path.islink(logfile):
            os.unlink(logfile)
    except Exception as e:
        print(f'Failed to delete {logfile}. Reason: {e}')

    arguments = ["scrapy",
                 "crawl",
                 spider,
                 "-s", f"LOG_FILE={logfile}",
                 "-s", f"staging_area={spider_sa_version}"]

    if version:
        arguments.append("-s")
        arguments.append(f"version={version}")

    if urls:
        arguments.append("-a")
        arguments.append(f"urls={urls}")

    return cmdline.execute(arguments)
