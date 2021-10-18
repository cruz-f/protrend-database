import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Source, Organism, RegulatoryFamily, Regulator, Operon, Gene, TFBS, Promoter, Effector, \
    RegulatoryInteraction
from protrend.transform.regulondb import EffectorTransformer
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.operon import OperonTransformer
from protrend.transform.regulondb.organism import OrganismTransformer
from protrend.transform.regulondb.promoter import PromoterTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.regulatory_family import RegulatoryFamilyTransformer
from protrend.transform.regulondb.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class SourceTransformer(RegulondbTransformer):
    name = 'regulondb'
    type = 'database'
    url = 'http://regulondb.ccg.unam.mx/'
    doi = '10.1093/nar/gky1077'
    authors = ['Alberto Santos-Zavaleta', 'Heladia Salgado', 'Socorro Gama-Castro', 'Mishael Sánchez-Pérez',
               'Laura Gómez-Romero', 'Daniela Ledezma-Tejeida', 'Jair Santiago García-Sotelo',
               'Kevin Alquicira-Hernández', 'Luis José Muñiz-Rascado', 'Pablo Peña-Loredo',
               'Cecilia Ishida-Gutiérrez', 'David A Velázquez-Ramírez', 'Víctor Del Moral-Chávez',
               'César Bonavides-Martínez', 'Carlos-Francisco Méndez-Cruz', 'James Galagan', 'Julio Collado-Vides']
    description = 'RegulonDB v 10.5: tackling challenges to unify classic and high throughput knowledge of gene regulation in E. coli K-12'

    default_node = Source
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self):
        db = dict(name=[self.name],
                        type=[self.type],
                        url=[self.url],
                        doi=[self.doi],
                        authors=[self.authors],
                        description=[self.description])

        df = pd.DataFrame(db, index=[0])

        self._stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(RegulondbConnector):
    default_from_node = Organism
    default_to_node = Source
    default_connect_stack = {'organism': 'integrated_organism.json', 'source': 'integrated_source.json'}

    def connect(self):
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = organism['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=['http://regulondb.ccg.unam.mx/'] * size,
                      external_identifier=['ECK12'] * size,
                      key=['organism'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryFamilyToSourceConnector(RegulondbConnector):
    default_from_node = RegulatoryFamily
    default_to_node = Source
    default_connect_stack = {'regulatory_family': 'integrated_regulatoryfamily.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rfam = read_from_stack(stack=self._connect_stack, file='regulatory_family',
                               default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = rfam['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for tf_id in rfam['transcription_factor_id']:
            if not is_null(tf_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={tf_id}&organism=ECK12&type=All')
                external_identifier.append(tf_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['transcription_factor_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToSourceConnector(RegulondbConnector):
    default_from_node = Regulator
    default_to_node = Source
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'source': 'integrated_source.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        keys = []
        for _, reg in regulator.iterrows():

            tf_id = reg.get('transcription_factor_id', None)
            srna_id = reg.get('srna_gene_id', None)
            sigma_id = reg.get('sigma_gene_id', None)

            if not is_null(tf_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={tf_id}&organism=ECK12&type=All')
                external_identifier.append(tf_id)
                keys.append('transcription_factor_id')
            elif not is_null(srna_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={srna_id}&organism=ECK12&type=All')
                external_identifier.append(srna_id)
                keys.append('srna_gene_id')
            elif not is_null(sigma_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={sigma_id}&organism=ECK12&type=All')
                external_identifier.append(sigma_id)
                keys.append('sigma_gene_id')
            else:
                url.append(None)
                external_identifier.append(None)
                keys.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=keys)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class OperonToSourceConnector(RegulondbConnector):
    default_from_node = Operon
    default_to_node = Source
    default_connect_stack = {'operon': 'integrated_operon.json', 'source': 'integrated_source.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for op_id in operon['operon_id']:
            if not is_null(op_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={op_id}&organism=ECK12&type=All')
                external_identifier.append(op_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['operon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class PromoterToSourceConnector(RegulondbConnector):
    default_from_node = Promoter
    default_to_node = Source
    default_connect_stack = {'promoter': 'integrated_promoter.json', 'source': 'integrated_source.json'}

    def connect(self):
        promoter = read_from_stack(stack=self._connect_stack, file='promoter',
                                   default_columns=PromoterTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = promoter['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for promoter_id in promoter['promoter_id']:
            if not is_null(promoter_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={promoter_id}&organism=ECK12&type=All')
                external_identifier.append(promoter_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['promoter_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class GeneToSourceConnector(RegulondbConnector):
    default_from_node = Gene
    default_to_node = Source
    default_connect_stack = {'gene': 'integrated_gene.json', 'source': 'integrated_source.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for gene_id in gene['gene_id']:
            if not is_null(gene_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={gene_id}&organism=ECK12&type=All')
                external_identifier.append(gene_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['gene_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class TFBSToSourceConnector(RegulondbConnector):
    default_from_node = TFBS
    default_to_node = Source
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'source': 'integrated_source.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for site_id in tfbs['site_id']:
            if not is_null(site_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={site_id}&organism=ECK12&type=All')
                external_identifier.append(site_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['site_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class EffectorToSourceConnector(RegulondbConnector):
    default_from_node = Effector
    default_to_node = Source
    default_connect_stack = {'effector': 'integrated_effector.json', 'source': 'integrated_source.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = effector['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for effector_id in effector['effector_id']:
            if not is_null(effector_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={effector_id}&organism=ECK12&type=All')
                external_identifier.append(effector_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['effector_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToSourceConnector(RegulondbConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://regulondb.ccg.unam.mx/search?term=ECK120011190&organism=ECK12&type=All
        url = []
        external_identifier = []
        for rin_id in rin['regulator_id']:
            if not is_null(rin_id):
                url.append(f'http://regulondb.ccg.unam.mx/search?term={rin_id}&organism=ECK12&type=All')
                external_identifier.append(rin_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['regulator_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
