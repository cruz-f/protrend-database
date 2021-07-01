from protrend.runners.spider_runner import run_spider

if __name__ == "__main__":
    # TODO: fix logfile
    run_spider(spider='regprecise',
               logfile='regprecise.log',
               user_name='neo4j',
               password='regprecise',
               ip='localhost',
               port='7687',
               db_name='neo4j',
               dbms=r'C:\Users\BiSBII\.Neo4jDesktop\relate-data\dbmss\dbms-9ebbc6bf-cbf0-456e-8fbe-3d57761bdcb8',
               import_folder=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src\protrend\extract\import\regprecise',
               version='0.0.0')