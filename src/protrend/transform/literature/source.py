import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Source, Organism, Regulator, Operon, Gene, RegulatoryInteraction, Effector
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.effector import EffectorTransformer
from protrend.transform.literature.gene import GeneTransformer
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.organism import OrganismTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.literature.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.processors import apply_processors, to_int_str
from protrend.utils import SetList


class SourceTransformer(LiteratureTransformer):
    name = ('bsub_faria_et_al_2017', 'ecol_fang_et_al_2017', 'mtub_turkarslan_et_al_2015', 'paer_vasquez_et_al_2011')
    type = ('literature', 'literature', 'literature', 'literature')
    url = ('https://www.frontiersin.org/articles/10.3389/fmicb.2016.00275',
           'https://www.pnas.org/content/114/38/10286',
           'https://www.nature.com/articles/sdata201510',
           'https://microbialinformaticsj.biomedcentral.com/articles/10.1186/2042-5783-1-3')
    doi = ('10.3389/fmicb.2016.00275',
           '10.1073/pnas.1702581114',
           '10.1038/sdata.2015.10',
           '10.1186/2042-5783-1-3')
    authors = (['José P. Faria', 'Ross Overbeek', 'Ronald C. Taylor', 'Neal Conrad', 'Veronika Vonstein',
                'Anne Goelzer', 'Vincent Fromion', 'Miguel Rocha', 'Isabel Rocha', 'Christopher S. Henry'],
               ['Xin Fang', 'Anand Sastry', 'Nathan Mih', 'Donghyuk Kim', 'Justin Tan', 'James T. Yurkovich',
                'Colton J. Lloyd', 'Ye Gao', 'Laurence Yang', 'Bernhard O. Palsson'],
               ['Serdar Turkarslan', 'Eliza J R Peterson', 'Tige R Rustad', 'Kyle J Minch', 'David J Reiss',
               'Robert Morrison', 'Shuyi Ma', 'Nathan D Price', 'David R Sherman', 'Nitin S Baliga'],
               ['Edgardo Galán-Vásquez', 'Beatriz Luna', 'Agustino Martínez-Antonio'])
    description = ('Reconstruction of the Regulatory Network for Bacillus subtilis and Reconciliation with Gene Expression Data',
                   'Global transcriptional regulatory network for Escherichia coli robustly connects gene expression to transcription factor activities',
                   'A comprehensive map of genome-wide gene regulation in Mycobacterium tuberculosis',
                   'The Regulatory Network of Pseudomonas aeruginosa')

    default_node = Source
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self):
        db = dict(name=self.name,
                  type=self.type,
                  url=self.url,
                  doi=self.doi,
                  authors=self.authors,
                  description=self.description)

        df = pd.DataFrame(db)

        self._stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(LiteratureConnector):
    default_from_node = Organism
    default_to_node = Source
    default_connect_stack = {'organism': 'integrated_organism.json', 'source': 'integrated_source.json'}

    def connect(self):
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        tax_to_source = {'224308': 'bsub_faria_et_al_2017',
                         '511145': 'ecol_fang_et_al_2017',
                         '83332': 'mtub_turkarslan_et_al_2015',
                         '208964': 'paer_vasquez_et_al_2011'}

        from_identifiers = []
        to_identifiers = []

        for organism_protrend_id, organism_id in zip(organism['protrend_id'], organism['ncbi_taxonomy']):
            from_identifiers.append(organism_protrend_id)

            to_id = tax_to_source[organism_id]
            to_identifiers.append(to_id)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToSourceConnector(LiteratureConnector):
    default_from_node = Regulator
    default_to_node = Source
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'source': 'integrated_source.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        df = pd.merge(regulator, source, left_on='source', right_on='name', suffixes=('_regulator', '_source'))

        from_identifiers = df['protrend_id_regulator'].tolist()
        to_identifiers = df['protrend_id_source'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToSourceConnector(LiteratureConnector):
    default_from_node = Operon
    default_to_node = Source
    default_connect_stack = {'operon': 'integrated_operon.json', 'source': 'integrated_source.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        df = pd.merge(operon, source, left_on='source', right_on='name', suffixes=('_operon', '_source'))

        from_identifiers = df['protrend_id_operon'].tolist()
        to_identifiers = df['protrend_id_source'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToSourceConnector(LiteratureConnector):
    default_from_node = Gene
    default_to_node = Source
    default_connect_stack = {'gene': 'integrated_gene.json', 'source': 'integrated_source.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        df = pd.merge(gene, source, left_on='source', right_on='name', suffixes=('_gene', '_source'))

        from_identifiers = df['protrend_id_gene'].tolist()
        to_identifiers = df['protrend_id_source'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EffectorToSourceConnector(LiteratureConnector):
    default_from_node = Effector
    default_to_node = Source
    default_connect_stack = {'effector': 'integrated_effector.json', 'source': 'integrated_source.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                               default_columns=EffectorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        df = pd.merge(effector, source, left_on='source', right_on='name', suffixes=('_effector', '_source'))

        from_identifiers = df['protrend_id_effector'].tolist()
        to_identifiers = df['protrend_id_source'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToSourceConnector(LiteratureConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        df = pd.merge(rin, source, left_on='source', right_on='name', suffixes=('_rin', '_source'))

        from_identifiers = df['protrend_id_rin'].tolist()
        to_identifiers = df['protrend_id_source'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)