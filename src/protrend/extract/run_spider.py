import os
import shutil

from scrapy import cmdline


def run_spider(spider: str,
               data_lake: str,
               version: str = None,
               urls: str = None):
    """


    Crawler settings
    :param spider:

    Spider pipeline
    :param data_lake:
    :param version:

    Spider settings
    :param urls:
    :return:
    """

    if not os.path.exists(data_lake):
        os.makedirs(data_lake)

    spider_dl = os.path.join(data_lake, spider)

    if not os.path.exists(spider_dl):
        os.makedirs(spider_dl)

    spider_dl_version = os.path.join(data_lake, spider, version)

    if not os.path.exists(spider_dl_version):
        os.makedirs(spider_dl_version)

    else:
        for filename in os.listdir(spider_dl_version):
            file_path = os.path.join(spider_dl_version, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')

    logfile = fr'{spider_dl_version}\{spider}.log'

    try:
        if os.path.isfile(logfile) or os.path.islink(logfile):
            os.unlink(logfile)
    except Exception as e:
        print(f'Failed to delete {logfile}. Reason: {e}')

    arguments = ["scrapy",
                 "crawl",
                 spider,
                 "-s", f"LOG_FILE={logfile}",
                 "-s", f"data_lake={spider_dl_version}"]

    if version:
        arguments.append("-s")
        arguments.append(f"version={version}")

    if urls:
        arguments.append("-a")
        arguments.append(f"urls={urls}")

    return cmdline.execute(arguments)
