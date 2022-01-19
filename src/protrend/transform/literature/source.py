from protrend.model import Source, Organism, Regulator, Gene, RegulatoryInteraction, Effector
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList
from protrend.utils.constants import LITERATURE
from protrend.utils.processors import to_int_str


class SourceTransformer(SourceMixIn, LiteratureTransformer,
                        source='literature',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['bsub_faria_et_al_2017', 'ecol_fang_et_al_2017', 'mtub_turkarslan_et_al_2015', 'paer_vasquez_et_al_2011']
    type = [LITERATURE, LITERATURE, LITERATURE, LITERATURE]
    url = ['https://www.frontiersin.org/articles/10.3389/fmicb.2016.00275',
           'https://www.pnas.org/content/114/38/10286',
           'https://www.nature.com/articles/sdata201510',
           'https://microbialinformaticsj.biomedcentral.com/articles/10.1186/2042-5783-1-3']
    doi = ['10.3389/fmicb.2016.00275',
           '10.1073/pnas.1702581114',
           '10.1038/sdata.2015.10',
           '10.1186/2042-5783-1-3']
    authors = [['José P. Faria', 'Ross Overbeek', 'Ronald C. Taylor', 'Neal Conrad', 'Veronika Vonstein',
                'Anne Goelzer', 'Vincent Fromion', 'Miguel Rocha', 'Isabel Rocha', 'Christopher S. Henry'],
               ['Xin Fang', 'Anand Sastry', 'Nathan Mih', 'Donghyuk Kim', 'Justin Tan', 'James T. Yurkovich',
                'Colton J. Lloyd', 'Ye Gao', 'Laurence Yang', 'Bernhard O. Palsson'],
               ['Serdar Turkarslan', 'Eliza J R Peterson', 'Tige R Rustad', 'Kyle J Minch', 'David J Reiss',
                'Robert Morrison', 'Shuyi Ma', 'Nathan D Price', 'David R Sherman', 'Nitin S Baliga'],
               ['Edgardo Galán-Vásquez', 'Beatriz Luna', 'Agustino Martínez-Antonio']]
    description = [
        'Reconstruction of the Regulatory Network for Bacillus subtilis and Reconciliation with Gene Expression Data',
        'Global transcriptional regulatory network for Escherichia coli robustly connects gene expression to '
        'transcription factor activities',
        'A comprehensive map of genome-wide gene regulation in Mycobacterium tuberculosis',
        'The Regulatory Network of Pseudomonas aeruginosa']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(LiteratureConnector,
                                source='literature',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='organism',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='name',
                                                     target_on='source',
                                                     source_processors={},
                                                     target_processors={'ncbi_taxonomy': [to_int_str]})

        taxa_to_source = {'224308': 'bsub_faria_et_al_2017',
                          '511145': 'ecol_fang_et_al_2017',
                          '83332': 'mtub_turkarslan_et_al_2015',
                          '208964': 'paer_vasquez_et_al_2011',
                          '1009714': 'paer_vasquez_et_al_2011',
                          '652611': 'paer_vasquez_et_al_2011',
                          '1081927': 'paer_vasquez_et_al_2011'}
        organism_source = target_df['ncbi_taxonomy'].map(taxa_to_source)
        target_df = target_df.assign(source=organism_source)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='name', target_on='source')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)


class SourceToRegulatorConnector(LiteratureConnector,
                                 source='literature',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self.create_connection(source='source', target='regulator', source_on='name', target_on='source')
        self.stack_connections(df)


class SourceToGeneConnector(LiteratureConnector,
                            source='literature',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self.create_connection(source='source', target='gene', source_on='name', target_on='source')
        self.stack_connections(df)


class SourceToEffectorConnector(LiteratureConnector,
                                source='literature',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Effector,
                                register=True):

    def connect(self):
        df = self.create_connection(source='source', target='effector', source_on='name', target_on='source')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(LiteratureConnector,
                                             source='literature',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self.create_connection(source='source', target='rin', source_on='name', target_on='source')
        self.stack_connections(df)
