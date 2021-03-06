from protrend.model import (Source, Organism, RegulatoryFamily, Regulator, Gene, TFBS, Effector,
                            RegulatoryInteraction)
from protrend.transform.mix_ins import SourceMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector
from protrend.utils import SetList, is_null
from protrend.utils.constants import DATABASE


class SourceTransformer(SourceMixIn, RegulonDBTransformer,
                        source='regulondb',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['regulondb']
    type = [DATABASE]
    url = ['http://regulondb.ccg.unam.mx/']
    doi = ['10.1093/nar/gky1077']
    authors = [['Alberto Santos-Zavaleta', 'Heladia Salgado', 'Socorro Gama-Castro', 'Mishael Sánchez-Pérez',
                'Laura Gómez-Romero', 'Daniela Ledezma-Tejeida', 'Jair Santiago García-Sotelo',
                'Kevin Alquicira-Hernández', 'Luis José Muñiz-Rascado', 'Pablo Peña-Loredo',
                'Cecilia Ishida-Gutiérrez', 'David A Velázquez-Ramírez', 'Víctor Del Moral-Chávez',
                'César Bonavides-Martínez', 'Carlos-Francisco Méndez-Cruz', 'James Galagan', 'Julio Collado-Vides']]
    description = ['RegulonDB v 10.5: tackling challenges to unify classic and high throughput knowledge of gene '
                   'regulation in E. coli K-12']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceConnector(RegulonDBConnector, register=False):

    def _connect(self, target: str, target_col: str, key: str):
        # noinspection DuplicatedCode
        source_df, target_df = self.transform_stacks(source='source',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        urls = []
        ids = []
        keys = []
        if target_col in target_df.columns:
            for target_id in target_df[target_col]:
                if not is_null(target_id):
                    urls.append(f'http://regulondb.ccg.unam.mx/search?term={target_id}&organism=ECK12&type=All')
                    ids.append(target_id)
                    keys.append(key)

                else:
                    urls.append(None)
                    ids.append(None)
                    keys.append(None)

        kwargs = dict(url=urls,
                      external_identifier=ids,
                      key=keys)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        return df


class SourceToOrganismConnector(RegulonDBConnector,
                                source='regulondb',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self.create_connection(source='source', target='organism')
        self.stack_connections(df)


class SourceToRegulatoryFamilyConnector(SourceConnector,
                                        source='regulondb',
                                        version='0.0.0',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):

    def connect(self):
        df = self._connect(target='rfam', target_col='transcription_factor_id', key='transcription_factor_id')
        self.stack_connections(df)


class SourceToRegulatorConnector(SourceConnector,
                                 source='regulondb',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self._connect(target='regulator', target_col='regulator_id', key='transcription_factor_id')
        self.stack_connections(df)


class SourceToGeneConnector(SourceConnector,
                            source='regulondb',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self._connect(target='gene', target_col='gene_id', key='gene_id')
        self.stack_connections(df)


class SourceToTFBSConnector(SourceConnector,
                            source='regulondb',
                            version='0.0.0',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):

    def connect(self):
        df = self._connect(target='tfbs', target_col='site_id', key='site_id')
        self.stack_connections(df)


class SourceToEffectorConnector(SourceConnector,
                                source='regulondb',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Effector,
                                register=True):

    def connect(self):
        df = self._connect(target='effector', target_col='effector_id', key='effector_id')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(SourceConnector,
                                             source='regulondb',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self._connect(target='rin', target_col='gene_id', key='gene_id')
        self.stack_connections(df)
