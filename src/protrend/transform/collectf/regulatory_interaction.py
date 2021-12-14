import pandas as pd

from protrend.io import read_json_frame, read_json_lines, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, Operon, Gene, TFBS
from protrend.transform import GeneTransformer
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_list, regulatory_effect_collectf, to_list_nan


class RegulatoryInteractionTransformer(CollectfTransformer,
                                       source='collectf',
                                       version='0.0.1',
                                       node=RegulatoryInteraction,
                                       order=70,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'integrated_tfbs.json',
                               'rin': 'TFBS.json'}
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_rin(self, rin: pd.DataFrame) -> pd.DataFrame:
        rin = self.select_columns(rin, 'tfbs_id', 'mode', 'regulon', 'gene')
        rin = apply_processors(rin, regulon=to_list_nan, gene=to_list_nan)
        rin = rin.explode('regulon')
        rin = rin.explode('gene')
        rin = rin.dropna(subset=['regulon', 'gene'])
        rin = rin.rename(columns={'mode': 'regulatory_effect'})
        return rin

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'tfbs_id', 'protrend_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'uniprot_accession', 'protrend_id', 'organism_protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator', 'organism_protrend_id': 'organism'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'locus_tag_old', 'protrend_id')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform(self) -> pd.DataFrame:
        rin = read_from_stack(stack=self.transform_stack, key='rin',
                              columns=TFBSTransformer.read_columns, reader=read_json_lines)

        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=TFBSTransformer.columns, reader=read_json_frame)

        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.read_columns, reader=read_json_frame)

        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.read_columns, reader=read_json_frame)

        rin = self.transform_rin(rin)

        # by regulator
        regulator = self.transform_regulator(regulator)
        rin_regulator = pd.merge(rin, regulator, left_on='regulon', right_on='uniprot_accession')

        # by gene
        gene = self.transform_gene(gene)
        rin_regulator_gene = pd.merge(rin_regulator, gene, left_on='gene', right_on='locus_tag_old')

        # by tfbs
        tfbs = self.transform_tfbs(tfbs)
        rin_regulator_gene_tfbs = pd.merge(rin_regulator_gene, tfbs, left_on='tfbs_id', right_on='tfbs_id')

        regulatory_interaction = apply_processors(rin_regulator_gene_tfbs,
                                                  regulatory_effect=regulatory_effect_collectf)
        regulatory_interaction = regulatory_interaction.assign(effector=None)

        regulatory_interaction = self.interaction_hash(regulatory_interaction)

        self.stack_transformed_nodes(regulatory_interaction)
        return regulatory_interaction


class RegulatoryInteractionToRegulatorConnector(CollectfConnector,
                                                source='collectf',
                                                version='0.0.1',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(CollectfConnector,
                                             source='collectf',
                                             version='0.0.1',
                                             from_node=RegulatoryInteraction,
                                             to_node=Operon,
                                             register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(CollectfConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_list)
        rin = rin.explode(column='genes')
        rin = rin.dropna(subset=['genes'])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['genes'].tolist()

        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(CollectfConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, tfbss=to_list)
        rin = rin.explode(column='tfbss')
        rin = rin.dropna(subset=['tfbss'])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['tfbss'].tolist()

        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToOperonConnector(CollectfConnector,
                                 source='collectf',
                                 version='0.0.1',
                                 from_node=Regulator,
                                 to_node=Operon,
                                 register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(CollectfConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_list)
        rin = rin.explode(column='genes')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['genes'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToTFBSConnector(CollectfConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, key='regulatory_interaction',
                              columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, tfbss=to_list)
        rin = rin.explode(column='tfbss')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['tfbss'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
