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
                       ])

    def transform(self) -> pd.DataFrame:
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'gene_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'regulator_effector'})

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

        # Genetic network
        # columns = 'regulator_id', 'regulator_name', 'regulated_id' 'regulated_name',
        # 'function_interaction', 'evidence', 'regulator_type', 'regulated_type'
        gen_net = self._build_genetic_network()

        tf_gen_net = pd.merge(gen_net, regulator, left_on='regulator_id', right_on='transcription_factor_id')
        sigma_gen_net = pd.merge(gen_net, regulator, left_on='regulator_id', right_on='sigma_id')

        gen_net = pd.concat([tf_gen_net, sigma_gen_net], axis=0)
        gen_net = pd.merge(gen_net, gene, left_on='regulated_id', right_on='gene_id')
        gen_net = gen_net.rename(columns={'regulated_id': 'gene_id', 'function_interaction': 'ri_function'})

        conformation = self._build_conformation()
        tf_conformation = pd.merge(conformation, regulator, on='transcription_factor_id')

        # Regulatory Interaction
        # columns = 'regulatory_interaction_id', 'conformation_id', 'promoter_id', 'site_id',
        # 'ri_function', 'ri_dist_first_gene'
        ri = self._build_regulatory_interaction()
        ri = pd.merge(ri, tf_conformation, on='conformation_id')
        ri = pd.merge(ri, gene, left_on='ri_dist_first_gene', right_on='gene_id')
        ri = ri.rename(columns={'transcription_factor_id': 'regulator_id'})

        # TF-Gene Interaction
        # columns = 'regulatory_interaction_id', 'conformation_id', 'object_id', 'site_id',
        # 'ri_function', 'ri_first_gene_id'
        interaction = self._build_tf_gene_interaction()
        interaction = pd.merge(interaction, tf_conformation, on='conformation_id')
        interaction = pd.merge(interaction, gene, left_on='object_id', right_on='gene_id')
        interaction = interaction.rename(columns={'transcription_factor_id': 'regulator_id'})

        # SRNA Interaction
        # columns = 'srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id', 'srna_function'
        srna = self._build_srna()
        srna = srna.rename(columns={'srna_gene_id': 'regulator_id', 'srna_gene_regulated_id': 'gene_id',
                                    'srna_function': 'ri_function'})

        regulator_gene = pd.concat([gen_net, ri, interaction, srna], axis=0)

        operon = self._build_operon()
        operon_by_gene = operon.explode(column='gene_id')
        regulator_operon = pd.merge(regulator_gene, operon_by_gene, on='gene_id')

        regulator_effector = self._build_regulator_effector()
        regulatory_interaction = pd.merge(regulator_operon, regulator_effector, how='left', on='regulator_id')
        regulatory_interaction = pd.merge(regulatory_interaction, effector, how='left', on='effector_id')

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

        rin = apply_processors(rin, effectors=to_list)
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
