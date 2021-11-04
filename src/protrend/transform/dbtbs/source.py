import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Source, Organism, RegulatoryFamily, Regulator, Operon, Gene, TFBS, \
    RegulatoryInteraction
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.dbtbs.organism import OrganismTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.regulatory_family import RegulatoryFamilyTransformer
from protrend.transform.dbtbs.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class SourceTransformer(DBTBSTransformer):
    name = 'dbtbs'
    type = 'database'
    url = 'https://dbtbs.hgc.jp/'
    doi = '10.1093/nar/gkm910'
    authors = ['Nicolas Sierro', 'Yuko Makita', 'Michiel de Hoon', 'Kenta Nakai']
    description = 'DBTBS: a database of transcriptional regulation in Bacillus subtilis containing upstream intergenic conservation information'

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


class OrganismToSourceConnector(DBTBSConnector):
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

        kwargs = dict(url=['https://dbtbs.hgc.jp/'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryFamilyToSourceConnector(DBTBSConnector):
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

        url = []
        external_identifier = []
        for rfam_name in rfam['name']:
            if not is_null(rfam_name) and rfam_name == 'Helix turn helix family':
                url.append('https://dbtbs.hgc.jp/tfactable.html#HTH')
                external_identifier.append('HTH')
            elif not is_null(rfam_name) and rfam_name == 'Sigma factors':
                url.append('https://dbtbs.hgc.jp/tfactable.html#sigma')
                external_identifier.append('sigma')
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['tfactable.html#'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToSourceConnector(DBTBSConnector):
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

        url = []
        external_identifier = []
        for _, reg in regulator.iterrows():

            reg_id = reg.get('name_dbtbs', None)
            reg_url = reg.get('url', None)

            if not is_null(reg_id) and not is_null(reg_url):
                url.append(reg_url)
                external_identifier.append(reg_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['tfac'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class OperonToSourceConnector(DBTBSConnector):
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

        url = []
        external_identifier = []
        for _, op in operon.iterrows():

            op_id = op.get('name', None)
            op_url = op.get('url', None)

            if not is_null(op_id) and not is_null(op_url):
                url.append(op_url)
                external_identifier.append(op_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['prom'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class GeneToSourceConnector(DBTBSConnector):
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

        url = []
        external_identifier = []
        for _, ge in gene.iterrows():

            ge_id = ge.get('name_dbtbs', None)
            ge_url = ge.get('url', None)

            if not is_null(ge_id) and not is_null(ge_url):
                url.append(ge_url)
                external_identifier.append(ge_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['prom'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class TFBSToSourceConnector(DBTBSConnector):
    default_from_node = TFBS
    default_to_node = Source
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'source': 'integrated_source.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.explode(column='operon')

        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        url = []
        external_identifier = []
        for _, bs in tfbs.iterrows():

            bs_id = bs.get('operon', None)
            bs_url = bs.get('url', None)

            if not is_null(bs_id) and not is_null(bs_url):
                url.append(bs_url)
                external_identifier.append(bs_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['prom'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToSourceConnector(DBTBSConnector):
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

        url = []
        external_identifier = []
        for _, ri in rin.iterrows():

            ri_id = ri.get('operon_name', None)
            ri_url = ri.get('url', None)

            if not is_null(ri_id) and not is_null(ri_url):
                url.append(ri_url)
                external_identifier.append(ri_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['prom'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
