import sys
sys.path.insert(0, r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src')


from protrend.runners.spider_runner import run_spider


if __name__ == "__main__":

    run_spider(spider='dbtbs',
               staging_area=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src\protrend\extract\staging_area',
               version='0.0.2')