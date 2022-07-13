import pandas as pd
from neo4j.exceptions import Neo4jError, DriverError

from protrend.transform import Transformer, Connector
from protrend.io import read_json_frame, read
from protrend.log import ProtrendLogger
from protrend.model import Gene, PromoterRegion, Organism
from protrend.transform.connector import CONNECTOR_STACK
from protrend.utils import Settings, apply_processors
from protrend.utils.processors import to_int_str


class PromoterRegionTransformer(Transformer,
                                source='promoter_region',
                                version='0.0.0',
                                node=PromoterRegion,
                                order=100,
                                register=True):

    @staticmethod
    def fetch_genes():
        try:
            return Gene.node_to_df()
        except (Neo4jError, DriverError) as e:
            ProtrendLogger.log.error(e)
            return pd.DataFrame()

    def transform(self) -> pd.DataFrame:
        genes = self.fetch_genes()
        genes = genes[['protrend_id', 'locus_tag']]

        promoters = []
        for genome_path in Settings.genomes_database.glob('*.json'):

            if not genome_path.exists():
                continue

            genome = read_json_frame(genome_path)
            genome = genome[['promoter_sequence', 'promoter_start', 'promoter_end', 'promoter_strand']]
            genome = genome.assign(ncbi_taxonomy=genome_path.stem)

            genome_promoters = pd.merge(genes, genome, on='locus_tag')
            promoters.append(genome_promoters)

        promoters = pd.concat(promoters, ignore_index=True)
        promoters = promoters.reset_index(drop=True)

        promoters = promoters.rename(columns={'protrend_id': 'gene',
                                              'promoter_sequence': 'sequence',
                                              'promoter_start': 'start',
                                              'promoter_end': 'stop',
                                              'promoter_strand': 'strand'})

        self.stack_transformed_nodes(promoters)
        return promoters


class PromoterRegionToOrganismConnector(Connector,
                                        source='promoter_region',
                                        version='0.0.0',
                                        from_node=PromoterRegion,
                                        to_node=Organism,
                                        register=True):

    @staticmethod
    def fetch_organisms():
        try:
            return Organism.node_to_df()
        except (Neo4jError, DriverError) as e:
            ProtrendLogger.log.error(e)
            return pd.DataFrame()

    def connect(self):
        file = CONNECTOR_STACK['promoter_region']
        source_df = read(source=self.source, version=self.version, file=file,
                         reader=read_json_frame, default=pd.DataFrame())
        source_df = apply_processors(source_df, ncbi_taxonomy=[to_int_str])
        source_df = source_df.rename(columns={'protrend_id': 'source_col'})

        target_df = self.fetch_organisms()
        target_df = apply_processors(target_df, ncbi_taxonomy=[to_int_str])
        target_df = target_df.rename(columns={'protrend_id': 'target_col'})

        df = pd.merge(source_df, target_df, on='ncbi_taxonomy')
        df = df.drop_duplicates(subset=['source_col', 'target_col'])

        source_ids = df['source_col'].to_list()
        target_ids = df['target_col'].to_list()

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)


class PromoterRegionToGeneConnector(Connector,
                                    source='promoter_region',
                                    version='0.0.0',
                                    from_node=PromoterRegion,
                                    to_node=Gene,
                                    register=True):

    def connect(self):
        df = self.create_connection(source='promoter_region', target='promoter_region',
                                    source_column='protrend_id', target_column='gene')
        self.stack_connections(df)
