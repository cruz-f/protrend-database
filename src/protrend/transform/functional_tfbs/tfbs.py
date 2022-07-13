import re

import pandas as pd
from neo4j.exceptions import Neo4jError, DriverError

from protrend.model import Regulator, Motif, TFBS, Organism
from protrend.transform.transformations import group_by
from protrend.binding_site_alignment import run_lasagna
from protrend.descriptors.consensus_sequence import consensus_sequence
from protrend.log import ProtrendLogger
from protrend.utils import Settings
from protrend.utils.constants import FORWARD, REVERSE
from protrend.utils.processors import take_last, to_list, to_list_nan
from protrend.transform import Transformer, Connector


NUCLEOTIDE_PATTERN = re.compile(r'[ATCG]+')


def get_aligned_sequence_start(aligned_sequence: str):
    if not aligned_sequence:
        return

    matches = list(re.finditer(NUCLEOTIDE_PATTERN, aligned_sequence))

    if not matches:
        return

    return matches[0].start()


def get_aligned_sequence_stop(aligned_sequence: str):
    if not aligned_sequence:
        return

    matches = list(re.finditer(NUCLEOTIDE_PATTERN, aligned_sequence))

    if not matches:
        return

    return matches[-1].end()


class TFBSTransformer(Transformer,
                      source='functional_tfbs',
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
                'strand': [],
                'ncbi_taxonomy': []}
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
            regulator_motif = motifs[motifs['regulator'] == regulator]
            regulator_motif = regulator_motif.reset_index(drop=True)

            headers = regulator_motif['tfbs'].to_list()
            sequences = regulator_motif['sequence'].to_list()

            aligned_motif = run_lasagna(headers, sequences, k=Settings.k)
            aligned_motif = pd.DataFrame(aligned_motif)

            aligned_motif = aligned_motif.assign(
                consensus_sequence=consensus_sequence(aligned_motif['aligned_sequence'].to_list())
            )

            motif = pd.merge(regulator_motif, aligned_motif, on='tfbs')
            motif = motif.reset_index(drop=True)

            aligned_motifs.append(motif)

        aligned_motifs = pd.concat(aligned_motifs, ignore_index=True)
        aligned_motifs = aligned_motifs.reset_index(drop=True)

        # remove the artifacts from lasagna
        aligned_motifs = aligned_motifs.drop(columns=['sequence'])
        aligned_motifs = aligned_motifs.rename(columns={'aligned_sequence': 'sequences'})

        aligned_motifs['strand'] = aligned_motifs['aligned_strand'].map({'+': FORWARD, '-': REVERSE})

        aligned_motifs['aligned_start'] = aligned_motifs['aligned_start'].map(get_aligned_sequence_start)
        aligned_motifs['aligned_stop'] = aligned_motifs['aligned_stop'].map(get_aligned_sequence_stop)
        aligned_motifs['start'] = aligned_motifs['start'] - aligned_motifs['aligned_start']
        aligned_motifs['stop'] = aligned_motifs['stop'] + aligned_motifs['aligned_stop']

        aligned_motifs = aligned_motifs.drop(columns=['aligned_start', 'aligned_stop', 'aligned_strand'])
        return aligned_motifs

    def transform(self) -> pd.DataFrame:
        motifs = self.fetch_motifs()
        aligned_motifs = self.align_motifs(motifs)

        aggregation = {'locus_tag': take_last,
                       'tfbs': to_list,
                       'organism': take_last,
                       'consensus_sequence': take_last,
                       'sequences': to_list,
                       'start': to_list,
                       'stop': to_list,
                       'strand': to_list}
        aligned_motifs = group_by(df=aligned_motifs,
                                  column='regulator',
                                  aggregation=aggregation, default=to_list)

        self.stack_transformed_nodes(aligned_motifs)
        return aligned_motifs


class MotifToTFBSConnector(Connector,
                           source='functional_tfbs',
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

        sequence = []
        start = []
        stop = []
        strand = []
        source_ids = []
        target_ids = []

        zipped = zip(source_df['protrend_id'], target_df['tfbs'], source_df['sequences'], source_df['start'],
                     source_df['stop'], source_df['strand'])
        for protrend_id, tfbs_ids, sequences, starts, stops, strands in zipped:
            for tfbs_id, sequence_, start_, stop_, strand_ in zip(tfbs_ids, sequences, starts, stops, strands):
                source_ids.append(protrend_id)
                target_ids.append(tfbs_id)
                sequence.append(sequence_)
                start.append(start_)
                stop.append(stop_)
                strand.append(strand_)

        kwargs = dict(sequence=sequence,
                      start=start,
                      stop=stop,
                      strand=strand)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class MotifToOrganismConnector(Connector,
                               source='functional_tfbs',
                               version='0.0.0',
                               from_node=Motif,
                               to_node=Organism,
                               register=True):

    def connect(self):
        df = self.create_connection(source='motif', target='motif',
                                    source_column='protrend_id', target_column='organism')
        self.stack_connections(df)


class MotifToRegulatorConnector(Connector,
                                source='functional_tfbs',
                                version='0.0.0',
                                from_node=Motif,
                                to_node=Regulator,
                                register=True):

    def connect(self):
        df = self.create_connection(source='motif', target='motif',
                                    source_column='protrend_id', target_column='regulator')
        self.stack_connections(df)
