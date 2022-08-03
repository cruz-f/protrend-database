import re

import pandas as pd
from neo4j.exceptions import Neo4jError, DriverError

from protrend.binding_site_alignment import run_lasagna
from protrend.descriptors.consensus_sequence import consensus_sequence
from protrend.log import ProtrendLogger
from protrend.model import Regulator, Motif, TFBS, Organism
from protrend.transform import Transformer, Connector
from protrend.transform.transformations import group_by, drop_empty_string, drop_duplicates
from protrend.utils import Settings
from protrend.utils.constants import FORWARD, REVERSE
from protrend.utils.processors import to_list, take_first

NUCLEOTIDE_PATTERN = re.compile(r'[ATCG]+')


def get_aligned_sequence_start(sequence_and_strand: tuple):
    aligned_sequence, strand = sequence_and_strand

    if not aligned_sequence:
        return

    matches = list(re.finditer(NUCLEOTIDE_PATTERN, aligned_sequence))

    if not matches:
        return

    if strand == FORWARD:
        return matches[0].start()

    elif strand == REVERSE:
        return matches[-1].end()

    else:
        return


def get_aligned_sequence_stop(sequence_and_strand: tuple):
    aligned_sequence, strand = sequence_and_strand

    if not aligned_sequence:
        return

    matches = list(re.finditer(NUCLEOTIDE_PATTERN, aligned_sequence))

    if not matches:
        return

    if strand == FORWARD:
        return matches[-1].end()

    elif strand == REVERSE:
        return matches[0].start()

    else:
        return


class MotifTransformer(Transformer,
                       source='motif',
                       version='0.0.0',
                       node=Motif,
                       order=100,
                       register=True):

    @staticmethod
    def fetch_motifs():
        try:
            regulators = Regulator.nodes.all()
        except (Neo4jError, DriverError) as e:
            regulators = []
            ProtrendLogger.log.error(e)

        data = {'regulator': [],
                'locus_tag': [],
                'tfbs': [],
                'organism': [],
                'sequence': [],
                'start': [],
                'stop': [],
                'strand': []}
        for regulator in regulators:
            sites = regulator.tfbs.all()
            if not sites:
                continue
            for site in sites:
                data['regulator'].append(regulator.protrend_id)
                data['locus_tag'].append(regulator.locus_tag)
                data['tfbs'].append(site.protrend_id)
                data['organism'].append(site.organism)
                data['sequence'].append(site.sequence)
                data['start'].append(site.start)
                data['stop'].append(site.stop)
                data['strand'].append(site.strand)

        df = pd.DataFrame(data)
        return df

    @staticmethod
    def align_motifs(motifs: pd.DataFrame) -> pd.DataFrame:
        """
        Aligns the motifs using lasagna.
        """
        aligned_motifs = []

        for regulator in motifs['regulator'].unique():
            regulator_motif = motifs[motifs['regulator'] == regulator].reset_index(drop=True)

            # it returns dict with 'identifiers', 'sequences', and 'strands'
            aligned_motif = run_lasagna(identifiers=regulator_motif['tfbs'].to_list(),
                                        sequences=regulator_motif['sequence'].to_list(),
                                        k=Settings.k)

            # we don't need the original 'sequence' anymore
            regulator_motif = regulator_motif.drop(columns=['sequence'])

            aligned_motif = pd.DataFrame({'tfbs': aligned_motif['identifiers'],
                                          'sequences': aligned_motif['sequences']})

            motif = pd.merge(regulator_motif, aligned_motif, on='tfbs')
            motif = motif.reset_index(drop=True)

            # creating the consensus sequence
            consensus = consensus_sequence(motif['sequences'])
            motif = motif.assign(consensus_sequence=consensus)

            aligned_motifs.append(motif)

        aligned_motifs = pd.concat(aligned_motifs, ignore_index=True)
        aligned_motifs = aligned_motifs.reset_index(drop=True)

        aligned_start = aligned_motifs[['sequences', 'strand']].apply(get_aligned_sequence_start, axis=1)
        aligned_motifs['start'] = aligned_motifs['start'] - aligned_start

        aligned_stop = aligned_motifs[['sequences', 'strand']].apply(get_aligned_sequence_stop, axis=1)
        aligned_motifs['stop'] = aligned_motifs['stop'] + aligned_stop

        return aligned_motifs

    def transform(self) -> pd.DataFrame:
        motifs = self.fetch_motifs()
        aligned_motifs = self.align_motifs(motifs)

        self.stack_transformed_nodes(aligned_motifs)
        return aligned_motifs

    def integrate(self, df: pd.DataFrame):
        """
        Integrates the motifs into the database.
        """
        # get the row duplicated aligned motifs dataframe
        aligned_motifs = df.copy()

        df = df[['locus_tag', 'regulator', 'tfbs', 'sequences', 'consensus_sequence']].copy()
        aggregation = {'locus_tag': take_first,
                       'tfbs': to_list,
                       'sequences': to_list,
                       'consensus_sequence': take_first}
        df = group_by(df=df,
                      column='regulator',
                      aggregation=aggregation, default=to_list)

        # take a db snapshot for the current node
        view = self.node_view()
        view = self.view_normalization(view)

        # ensure uniqueness
        df, factorized_cols = self.factors_normalization(df=df)

        # assign the integration and load columns
        df = df.assign(protrend_id=None, load=None, what='nodes')

        # try to integrate the new nodes by the node factors
        for factor in factorized_cols:
            # It should not be the case, but there may be some duplicated records in the database and empty values.
            # So we must prevent this
            mapper = view.dropna(subset=[factor])
            mapper = drop_empty_string(mapper, factor)
            mapper = drop_duplicates(mapper, subset=[factor])
            mapper = mapper.set_index(factor)
            mapper = mapper['protrend_id']

            to_fill = df[factor].map(mapper)
            filled = df['protrend_id'].fillna(to_fill)
            df = df.assign(protrend_id=filled)

        # those that mismatch the integration by the node factors (that is, the protrend id is missing) will be assigned
        # for creation
        mask = df['protrend_id'].isna()
        df.loc[mask, 'load'] = 'create'
        df.loc[~mask, 'load'] = 'update'

        # inferring the batch size, creating the new ids and assign to the df
        nodes = view['protrend_id'].to_list()
        last_node_idx = self.last_node_index(nodes)
        batch_size = mask.sum()

        batch_ids = self.protrend_identifiers_batch(last_node_idx=last_node_idx, size=batch_size)
        df.loc[mask, 'protrend_id'] = batch_ids

        # get the protrend ids in the integrated dataframe
        integrated_df = df[['protrend_id', 'regulator']]
        aligned_motifs = aligned_motifs.merge(integrated_df, on='regulator')

        self.stack_integrated_nodes(aligned_motifs)
        self.stack_nodes(df)


class MotifToTFBSConnector(Connector,
                           source='motif',
                           version='0.0.0',
                           from_node=Motif,
                           to_node=TFBS,
                           register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='motif',
                                                     target='motif',
                                                     source_column='protrend_id',
                                                     target_column='tfbs',
                                                     source_processors={},
                                                     target_processors={})

        source_df = source_df.drop_duplicates(subset=['source_col', 'tfbs'])

        source_ids = source_df['source_col'].to_list()
        target_ids = source_df['tfbs'].to_list()

        kwargs = dict(sequence=source_df['sequences'].to_list(),
                      start=source_df['start'].to_list(),
                      stop=source_df['stop'].to_list(),
                      strand=source_df['strand'].to_list())

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class MotifToOrganismConnector(Connector,
                               source='motif',
                               version='0.0.0',
                               from_node=Motif,
                               to_node=Organism,
                               register=True):

    def connect(self):
        df = self.create_connection(source='motif', target='motif',
                                    source_column='protrend_id', target_column='organism')
        self.stack_connections(df)


class MotifToRegulatorConnector(Connector,
                                source='motif',
                                version='0.0.0',
                                from_node=Motif,
                                to_node=Regulator,
                                register=True):

    def connect(self):
        df = self.create_connection(source='motif', target='motif',
                                    source_column='protrend_id', target_column='regulator')
        self.stack_connections(df)
