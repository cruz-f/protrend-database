import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction, Regulator, TFBS, Gene, Operon, Effector
from protrend.transform.processors import (apply_processors, to_list_nan, regulatory_effect_regulondb, to_list)
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.effector import EffectorTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.operon import OperonTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.utils import SetList


class RegulatoryInteractionTransformer(RegulondbTransformer):
    default_node = RegulatoryInteraction
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json',
                               'effector': 'integrated_effector.json',
                               'gene': 'integrated_gene.json'}
    default_order = 70
    columns = SetList(['regulator', 'operon', 'genes', 'tfbss', 'regulator_effector', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'gene_id', 'effector_id', 'transcription_factor_id', 'srna_gene_id', 'sigma_id', 'operon_id',
                       'regulator_id', 'evidence'])

    def transform(self) -> pd.DataFrame:
        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'regulator_effector'})

        # columns = ['regulator_id', 'effector_id', 'regulator_effector']
        regulator_effector = self._build_regulator_effector()
        effector = pd.merge(regulator_effector, effector, on='effector_id')

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator,
                                        'transcription_factor_id',
                                        'srna_gene_id',
                                        'sigma_id',
                                        'protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'genes', 'tfbss', 'operon_id', 'protrend_id')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list_nan, tfbss=to_list_nan)
        operon_by_gene = operon.explode(column='genes')

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'gene_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        # columns = ['genes', 'tfbss', 'operon_id', 'operon', 'gene_id', 'gene_protrend_id']
        operon_by_gene = pd.merge(operon_by_gene, gene, left_on='genes', right_on='gene_protrend_id')

        # columns = ['regulator_id', 'gene_id', 'site_id', 'ri_function', 'evidence',
        # 'genes', 'tfbss', 'operon_id', 'operon', 'gene_protrend_id']
        meta_regulatory_interaction = self._build_meta_regulatory_interaction()
        meta_regulatory_interaction = pd.merge(meta_regulatory_interaction, operon_by_gene, on='gene_id')

        regulator_tf = regulator.dropna(subset=['transcription_factor_id'])
        regulatory_interaction_tf = pd.merge(meta_regulatory_interaction, regulator_tf,
                                             left_on='regulator_id', right_on='transcription_factor_id')

        regulator_srna = regulator.dropna(subset=['srna_gene_id'])
        regulatory_interaction_srna = pd.merge(meta_regulatory_interaction, regulator_srna,
                                               left_on='regulator_id', right_on='srna_gene_id')

        regulator_sigma = regulator.dropna(subset=['sigma_id'])
        regulatory_interaction_sigma = pd.merge(meta_regulatory_interaction, regulator_sigma,
                                                left_on='regulator_id', right_on='sigma_id')

        # columns = ['regulator_id', 'gene_id', 'site_id', 'ri_function', 'evidence',
        # 'genes', 'tfbss', 'operon_id', 'operon', 'gene_protrend_id', 'transcription_factor_id', 'srna_gene_id',
        # 'sigma_id', 'regulator']
        regulatory_interaction = pd.concat([regulatory_interaction_tf,
                                            regulatory_interaction_srna,
                                            regulatory_interaction_sigma], axis=0)

        # add effectors
        regulatory_interaction = pd.merge(regulatory_interaction, effector, how='left', on='regulator_id')

        # add genes and tfbs
        regulatory_interaction = regulatory_interaction.drop(columns=['genes', 'tfbss', 'operon_id'])
        regulatory_interaction = pd.merge(regulatory_interaction, operon, on='operon')

        regulatory_interaction = regulatory_interaction.rename(columns={'ri_function': 'regulatory_effect'})
        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_regulondb)

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToEffectorConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Effector
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, regulator_effector=to_list)
        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator_effector'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Regulator
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Operon
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Gene
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

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


class RegulatoryInteractionToTFBSConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = TFBS
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

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


class RegulatorToEffectorConnector(RegulondbConnector):
    default_from_node = Regulator
    default_to_node = Effector
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = rin.drop_duplicates(subset=['regulator', 'regulator_effector'])

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['regulator_effector'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToOperonConnector(RegulondbConnector):
    default_from_node = Regulator
    default_to_node = Operon
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(RegulondbConnector):
    default_from_node = Regulator
    default_to_node = Gene
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_list)
        rin = rin.explode(column='genes')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['genes'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToTFBSConnector(RegulondbConnector):
    default_from_node = Regulator
    default_to_node = TFBS
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, tfbss=to_list)
        rin = rin.explode(column='tfbss')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['tfbss'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
