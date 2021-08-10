import re
from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.transform.processors import apply_processors, remove_ellipsis, \
    upper_case, tfbs_left_position, operon_left_position, operon_strand, tfbs_right_position
from protrend.transform.regprecise.settings import TFBSSettings
from protrend.transform.transformer import Transformer

regprecise_tfbs_pattern = re.compile(r'-\([0-9]+\)-')

# 'protrend_id',
#                                        'position', 'score', 'sequence',
#                                        'tfbs_id', 'url', 'regulon',
#                                        'operon', 'gene', 'tfbs_id_old',
#                                        'position_left', 'position_right'
class TFBSTransformer(Transformer):

    def __init__(self, settings: TFBSSettings = None):

        if not settings:
            settings = TFBSSettings()

        super().__init__(settings)

    def _read_tfbs(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('tfbs')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['position', 'score', 'sequence',
                                       'tfbs_id', 'url', 'regulon',
                                       'operon', 'gene'])

        return df

    def _read_gene(self) -> pd.DataFrame:

        file_path = self._transform_stack.get('gene')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=['protrend_id',
                                       'url', 'regulon', 'operon', 'tfbs',
                                       'organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
                                       'regulator_protrend_id', 'regulon_id',
                                       'locus_tag', 'name',
                                       'synonyms', 'function',
                                       'description'
                                       'ncbi_gene', 'ncbi_protein',
                                       'genbank_accession', 'refseq_accession', 'uniprot_accession',
                                       'sequence',
                                       'strand',
                                       'position_left',
                                       'position_right',
                                       'annotation_score'
                                       'locus_tag_regprecise'])

        df = df.dropna(subset=['protrend_id'])
        df = df.dropna(subset=['locus_tag_regprecise'])
        df = df.drop_duplicates(subset=['locus_tag_regprecise'])

        df = df[['strand', 'position_left', 'locus_tag_regprecise']]

        df = df.set_index(df['locus_tag_regprecise'])

        return df

    @staticmethod
    def _reduce_sequence(position, sequence):

        all_subgroups: List[Union[re.Match, None]] = list(re.finditer(regprecise_tfbs_pattern, sequence))

        all_subgroups: List[Union[re.Match, None]] = all_subgroups + [None]

        sequences = []
        last_pos = position
        start = 0

        for subgroup in all_subgroups:

            if subgroup is None:
                seq = sequence[start:len(sequence)]
                pos = last_pos
                sequences.append((seq, pos))

            else:
                group = subgroup.group().replace('-(', '').replace(')-', '')
                end = subgroup.start()

                seq = sequence[start:end]
                pos = last_pos
                sequences.append((seq, pos))

                start = subgroup.end()

                length = int(group)
                last_pos = last_pos + len(seq) + length

        return sequences

    @staticmethod
    def _new_tfbs_identifier(tfbs_id: str, pos: int) -> str:

        _, *genes = tfbs_id.split('_')

        tfbs = [str(pos)] + list(genes)

        return '_'.join(tfbs)

    def _normalize_sequence(self, df: pd.DataFrame) -> pd.DataFrame:

        res = defaultdict(list)

        for _, row in df.iterrows():

            position = row['position']
            score = row['score']
            sequence = row['sequence']
            tfbs_id = row['tfbs_id']
            url = row['url']
            regulon = row['regulon']
            operon = row['operon']
            gene = row['gene']

            new_seqs = self._reduce_sequence(position=position, sequence=sequence)

            for seq, pos in new_seqs:
                new_tfbs_id = self._new_tfbs_identifier(tfbs_id, pos)

                res['position'].append(pos)
                res['score'].append(score)
                res['sequence'].append(seq)
                res['tfbs_id_old'].append(tfbs_id)
                res['tfbs_id'].append(new_tfbs_id)
                res['url'].append(url)
                res['regulon'].append(regulon)
                res['operon'].append(operon)
                res['gene'].append(gene)

        return pd.DataFrame(res)

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:

        tfbs = tfbs.drop_duplicates(subset=['tfbs_id'])
        tfbs = tfbs.dropna(subset=['sequence'])
        tfbs = tfbs.reset_index(drop=True)

        apply_processors(remove_ellipsis,
                         df=tfbs,
                         col='sequence')

        apply_processors(upper_case,
                         df=tfbs,
                         col='sequence')

        tfbs = self._normalize_sequence(tfbs)
        tfbs = tfbs.explode('regulon')
        tfbs = tfbs.drop_duplicates(subset=['position', 'sequence', 'tfbs_id', 'regulon'])

        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        strands = []
        positions_left = []
        positions_right = []

        for position, seq, genes in zip(tfbs['position'], tfbs['sequence'], tfbs['gene']):

            op_strand = None

            for ge in genes:
                ge_strand = gene.loc[ge, 'strand']
                op_strand = operon_strand(previous_strand=op_strand,
                                          current_strand=ge_strand)

            operon_left = None

            for ge in genes:
                ge_left = gene.loc[ge, 'position_left']

                operon_left = operon_left_position(strand=op_strand,
                                                   previous_left=operon_left,
                                                   current_left=ge_left)

            tfbs_left = tfbs_left_position(strand=op_strand,
                                           gene_position=operon_left,
                                           gene_relative_position=position)

            tfbs_right = tfbs_right_position(strand=op_strand,
                                             gene_position=operon_left,
                                             gene_relative_position=position,
                                             tfbs_length=len(seq))

            strands.append(op_strand)
            positions_left.append(tfbs_left)
            positions_right.append(tfbs_right)

        tfbs['position_left'] = positions_left
        tfbs['position_right'] = positions_right

        return tfbs

    def transform(self):
        tfbs = self._read_tfbs()
        tfbs = self._transform_tfbs(tfbs)

        gene = self._read_gene()

        df = self._tfbs_coordinates(tfbs, gene)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def connect(self):

        pass
