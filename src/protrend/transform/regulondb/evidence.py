from typing import Union

import pandas as pd

from protrend.io import read_from_stack, read_txt, read_json_frame
from protrend.model import Evidence, Promoter, Regulator, TFBS, Operon, Gene, RegulatoryInteraction
from protrend.utils.processors import apply_processors, rstrip, lstrip, split_semi_colon, to_list_nan
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.operon import OperonTransformer
from protrend.transform.regulondb.promoter import PromoterTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.transform.regulondb.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils import SetList


class EvidenceTransformer(RegulondbTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Evidence,
                          order=100,
                          register=True):
    default_transform_stack = {'evidence': 'evidence.txt'}
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                       'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                       'example'])
    read_columns = SetList(['evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                            'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                            'example'])

    def _transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=evidence, subset=['evidence_id', 'evidence_name'])

        def remove_evidence_note(item: str) -> Union[None, str]:
            item = item.replace('EVIDENCE_NOTE', '')
            if not item:
                return
            return item

        df = apply_processors(df,
                              evidence_id=[rstrip, lstrip], evidence_name=[rstrip, lstrip],
                              evidence_note=[rstrip, lstrip, remove_evidence_note])
        df = df.dropna(subset=['evidence_id', 'evidence_name'])
        df['name'] = df['evidence_name']
        df['description'] = df['evidence_note']

        return df

    def transform(self):
        evidence = read_from_stack(stack=self.transform_stack, file='evidence', default_columns=self.read_columns,
                                   reader=read_txt, skiprows=38, names=self.read_columns)
        evidence = self._transform_evidence(evidence)

        self._stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatorConnector(RegulondbConnector,
                                   source='regulondb',
                                   version='0.0.0',
                                   from_node=Evidence,
                                   to_node=Regulator,
                                   register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'regulator': 'integrated_regulator.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'evidence_id']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator[['protrend_id', 'transcription_factor_id', 'sigma_id', 'srna_id']]
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        regulator_tf = regulator.dropna(subset=['transcription_factor_id'])
        obj_tf = pd.merge(obj_ev_pub, regulator_tf, left_on='object_id', right_on='transcription_factor_id')

        regulator_sigma = regulator.dropna(subset=['sigma_id'])
        obj_sigma = pd.merge(obj_ev_pub, regulator_sigma, left_on='object_id', right_on='sigma_id')

        regulator_srna = regulator.dropna(subset=['srna_id'])
        obj_srna = pd.merge(obj_ev_pub, regulator_srna, left_on='object_id', right_on='srna_id')

        obj_reg = pd.concat([obj_tf, obj_sigma, obj_srna])

        reg_ev = pd.merge(obj_reg, evidence, on='evidence_id')
        reg_ev = reg_ev.drop_duplicates(subset=['evidence_protrend_id', 'regulator_protrend_id'])

        from_identifiers = reg_ev['evidence_protrend_id'].tolist()
        to_identifiers = reg_ev['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToTFBSConnector(RegulondbConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=Evidence,
                              to_node=TFBS,
                              register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'tfbs': 'integrated_tfbs.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'evidence_id']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs[['protrend_id', 'site_id']]
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_site = pd.merge(obj_ev_pub, tfbs, left_on='object_id', right_on='site_id')

        tfbs_ev = pd.merge(obj_site, evidence, on='evidence_id')
        tfbs_ev = tfbs_ev.drop_duplicates(subset=['evidence_protrend_id', 'tfbs_protrend_id'])

        from_identifiers = tfbs_ev['evidence_protrend_id'].tolist()
        to_identifiers = tfbs_ev['tfbs_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToPromoterConnector(RegulondbConnector,
                                  source='regulondb',
                                  version='0.0.0',
                                  from_node=Evidence,
                                  to_node=Promoter,
                                  register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'promoter': 'integrated_promoter.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'evidence_id']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        promoter = read_from_stack(stack=self._connect_stack, file='promoter',
                                   default_columns=PromoterTransformer.columns, reader=read_json_frame)
        promoter = promoter[['protrend_id', 'promoter_id']]
        promoter = promoter.rename(columns={'protrend_id': 'promoter_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_promoter = pd.merge(obj_ev_pub, promoter, left_on='object_id', right_on='promoter_id')

        promoter_ev = pd.merge(obj_promoter, evidence, on='evidence_id')
        promoter_ev = promoter_ev.drop_duplicates(subset=['evidence_protrend_id', 'promoter_protrend_id'])

        from_identifiers = promoter_ev['evidence_protrend_id'].tolist()
        to_identifiers = promoter_ev['promoter_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToOperonConnector(RegulondbConnector,
                                source='regulondb',
                                version='0.0.0',
                                from_node=Evidence,
                                to_node=Operon,
                                register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'operon': 'integrated_operon.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'evidence_id']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = operon[['protrend_id', 'operon_id']]
        operon = operon.rename(columns={'protrend_id': 'operon_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_operon = pd.merge(obj_ev_pub, operon, left_on='object_id', right_on='operon_id')

        operon_ev = pd.merge(obj_operon, evidence, on='evidence_id')
        operon_ev = operon_ev.drop_duplicates(subset=['evidence_protrend_id', 'operon_protrend_id'])

        from_identifiers = operon_ev['evidence_protrend_id'].tolist()
        to_identifiers = operon_ev['operon_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToGeneConnector(RegulondbConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=Evidence,
                              to_node=Gene,
                              register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'gene': 'integrated_gene.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'evidence_id']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = gene[['protrend_id', 'gene_id']]
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_gene = pd.merge(obj_ev_pub, gene, left_on='object_id', right_on='gene_id')

        gene_ev = pd.merge(obj_gene, evidence, on='evidence_id')
        gene_ev = gene_ev.drop_duplicates(subset=['evidence_protrend_id', 'gene_protrend_id'])

        from_identifiers = gene_ev['evidence_protrend_id'].tolist()
        to_identifiers = gene_ev['gene_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToRegulatoryInteractionConnector(RegulondbConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=Evidence,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence[['protrend_id', 'name']]
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = rin[['protrend_id', 'evidence']]
        rin = rin.rename(columns={'protrend_id': 'regulatory_interaction_protrend_id'})
        rin = apply_processors(rin, evidence=[rstrip, lstrip, split_semi_colon, to_list_nan])
        rin = rin.explode(column='evidence')

        df = pd.merge(evidence, rin, left_on='name', right_on='evidence')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'regulatory_interaction_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['regulatory_interaction_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
