import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction
from protrend.transform.processors import (apply_processors, to_list, to_list_nan, regulatory_interaction_hash,
                                           regulatory_effect_regulondb)
from protrend.transform.regulondb.base import RegulondbTransformer
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
    columns = SetList(['operon', 'regulon', 'operon_id', 'regulatory_effect', 'regulator',
                       'organism_protrend_id', 'genes', 'tfbss',
                       'regulatory_interaction_hash', 'protrend_id'])

    def transform(self) -> pd.DataFrame:
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'gene_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effectors'})

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
        gen_net = self._build_genetic_network()

        tf_gen_net = pd.merge(gen_net, regulator, left_on='regulator_id', right_on='transcription_factor_id')
        sigma_gen_net = pd.merge(gen_net, regulator, left_on='regulator_id', right_on='sigma_id')

        gen_net = pd.concat([tf_gen_net, sigma_gen_net], axis=0)
        gen_net = pd.merge(gen_net, gene, left_on='regulated_id', right_on='gene_id')
        gen_net = gen_net.rename(columns={'regulated_id': 'gene_id'})

        conformation = self._build_conformation()
        tf_conformation = pd.merge(conformation, regulator, on='transcription_factor_id')

        # Regulatory Interaction
        ri = self._build_regulatory_interaction()
        ri = pd.merge(ri, tf_conformation, on='conformation_id')
        ri = pd.merge(ri, gene, left_on='ri_dist_first_gene', right_on='gene_id')
        ri = ri.rename(columns={'transcription_factor_id': 'regulator_id'})

        # TF-Gene Interaction
        interaction = self._build_tf_gene_interaction()
        interaction = pd.merge(interaction, tf_conformation, on='conformation_id')
        interaction = pd.merge(interaction, gene, left_on='object_id', right_on='gene_id')
        interaction = interaction.rename(columns={'transcription_factor_id': 'regulator_id'})

        # SRNA Interaction
        srna = self._build_srna()
        srna = srna.rename(columns={'srna_gene_id': 'regulator_id', 'srna_gene_regulated_id': 'gene_id'})

        regulator_gene = pd.concat([gen_net, ri, interaction, srna], axis=0)

        operon = self._build_operon()
        operon_by_gene = operon.explode(column='gene_id')
        regulator_operon = pd.merge(regulator_gene, operon_by_gene, on='gene_id')

        regulator_effector = self._build_regulator_effector()
        regulatory_interaction = pd.merge(regulator_operon, regulator_effector, how='left', on='regulator_id')

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  function_interaction=regulatory_effect_regulondb)
        regulatory_interaction = regulatory_interaction.rename(columns={'function_interaction': 'regulatory_effect'})

        # filter by regulatory effect + effector + regulator + operon
        df2 = apply_processors(regulatory_interaction,
                               regulatory_effect=to_list, effector=to_list, regulator=to_list, operon=to_list)

        _hash = df2['effectors'] + df2['regulator'] + df2['operon'] + df2['regulatory_effect']
        regulatory_interaction['regulatory_interaction_hash'] = _hash

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_interaction_hash=regulatory_interaction_hash)
        regulatory_interaction = regulatory_interaction.dropna(subset=['regulatory_interaction_hash'])

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction
