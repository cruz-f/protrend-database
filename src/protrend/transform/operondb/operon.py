import pandas as pd

from protrend.io import read_txt
from protrend.io.utils import read_gene, read
from protrend.model import Operon, Gene, Organism
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.transform.operondb.gene import GeneTransformer
from protrend.transform.transformations import select_columns, group_by, drop_empty_string
from protrend.utils import SetList
from protrend.utils.constants import FORWARD, REVERSE
from protrend.utils.processors import to_set_list, take_last, strand_mode, start_forward, start_reverse, to_list_nan, \
    apply_processors


class OperonTransformer(OperonDBTransformer,
                        source='operondb',
                        version='0.0.0',
                        node=Operon,
                        order=90,
                        register=True):
    columns = SetList(['protrend_id', 'operon_db_id', 'name', 'function', 'genes', 'strand', 'start', 'stop',
                       'organism', 'source'])

    @staticmethod
    def transform_gene(operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'strand', 'start', 'stop',
                              'organism', 'operon_db_id')
        gene = apply_processors(gene, operon_db_id=to_list_nan)
        gene = gene.explode('operon_db_id')
        gene = gene.rename(columns={'protrend_id': 'genes'})

        operon = select_columns(operon, 'operon_db_id', 'name', 'definition', 'source')
        operon = operon.rename(columns={'definition': 'function'})

        operon = pd.merge(operon, gene, on='operon_db_id')

        aggregation = {'genes': to_set_list, 'strand': to_set_list, 'start': to_set_list, 'stop': to_set_list,
                       'name': take_last, 'function': take_last, 'source': take_last, 'organism': take_last}
        operon = group_by(operon, column='operon_db_id', aggregation=aggregation)

        operon = operon.dropna(subset=['operon_db_id'])
        operon = drop_empty_string(operon, 'operon_db_id')
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
        conserved_columns = ['coid', 'org', 'name', 'op', 'definition', 'source', 'mbgd']
        known_columns = ['koid', 'org', 'name', 'op', 'definition', 'source']
        conserved = read(source=self.source, version=self.version, file='conserved_operon.txt', reader=read_txt,
                         default=pd.DataFrame(columns=conserved_columns))
        known = read(source=self.source, version=self.version, file='known_operon.txt', reader=read_txt,
                     default=pd.DataFrame(columns=known_columns))

        operon = self.transform_operon(conserved=conserved, known=known)

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

        operon = self.transform_gene(operon, gene)
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
                                                     target='gene',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='operon_db_id',
                                                     target_on='operon_db_id',
                                                     source_processors={},
                                                     target_processors={'operon_db_id': [to_list_nan]})
        target_df = target_df.explode('operon_db_id')

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
