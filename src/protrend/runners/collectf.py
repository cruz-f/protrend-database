from protrend.runners.spider_runner import run_spider

if __name__ == "__main__":
    # TODO: fix logfile

    run_spider(spider='collectf',
               logfile='collectf.log',
               user_name='neo4j',
               password='collectf',
               ip='localhost',
               port='7687',
               db_name='neo4j',
               dbms=r'C:\Users\BiSBII\.Neo4jDesktop\relate-data\dbmss\dbms-1a131eeb-1263-417b-8535-cdb9a0d10811',
               import_folder=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src\protrend\extract\import\collectf',
               version='0.0.0')