import pandas as pd

from protrend.io.utils import read_gene
from protrend.model import Operon, Gene, Organism
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.transform.operondb.gene import GeneTransformer
from protrend.transform.transformations import select_columns, group_by
from protrend.utils import SetList
from protrend.utils.constants import FORWARD, REVERSE
from protrend.utils.processors import to_set_list, take_last, strand_mode, start_forward, start_reverse, to_list_nan


class OperonTransformer(OperonDBTransformer,
                        source='operondb',
                        version='0.0.0',
                        node=Operon,
                        order=90,
                        register=True):
    columns = SetList(['protrend_id', 'operon_db_id', 'name', 'function', 'genes', 'strand', 'start', 'stop',
                       'organism', 'pubmed'])

    @staticmethod
    def transform_gene(gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene,
                              'protrend_id', 'strand', 'start', 'stop',
                              'organism',
                              'operon_db_id', 'operon_name', 'operon_function', 'pubmed')
        gene = gene.rename(columns={'protrend_id': 'genes', 'operon_name': 'name', 'operon_function': 'function'})

        aggregation = {'genes': to_set_list, 'strand': to_set_list, 'start': to_set_list, 'stop': to_set_list,
                       'name': take_last, 'function': take_last, 'pubmed': take_last, 'organism': take_last}
        operon = group_by(gene, column='operon_db_id', aggregation=aggregation)
        return operon

    @staticmethod
    def operon_coordinates(operon: pd.DataFrame) -> pd.DataFrame:
        strands = []
        starts = []
        stops = []
        for strand, start, stop in zip(operon['strand'], operon['start'], operon['stop']):

            strand = strand_mode(strand)

            if strand == FORWARD:
                start = start_forward(start)
                stop = start_reverse(stop)

            elif strand == REVERSE:
                start = start_reverse(start)
                stop = start_forward(stop)

            else:
                start = None
                stop = None

            strands.append(strand)
            starts.append(start)
            stops.append(stop)

        operon = operon.assign(strand=strands, start=starts, stop=stops)
        return operon

    def transform(self):
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

        operon = self.transform_gene(gene)
        operon = self.operon_coordinates(operon)

        self.stack_transformed_nodes(operon)
        return operon


class OperonToGeneConnector(OperonDBConnector,
                            source='operondb',
                            version='0.0.0',
                            from_node=Operon,
                            to_node=Gene,
                            register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='operon',
                                                     target='operon',
                                                     source_column='protrend_id',
                                                     target_column='genes',
                                                     source_on='operon_db_id',
                                                     target_on='operon_db_id',
                                                     source_processors={},
                                                     target_processors={'genes': [to_list_nan]})
        target_df = target_df.explode('genes')

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='operon_db_id', target_on='operon_db_id')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)


class OperonToOrganismConnector(OperonDBConnector,
                                source='operondb',
                                version='0.0.0',
                                from_node=Operon,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self.create_connection(source='operon', target='operon',
                                    target_column='organism')
        self.stack_connections(df)
