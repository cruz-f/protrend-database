import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Source, Organism, RegulatoryFamily, Regulator, Operon, Gene, TFBS, \
    RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.operon import OperonTransformer
from protrend.transform.collectf.organism import OrganismTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.regulatory_family import RegulatoryFamilyTransformer
from protrend.transform.collectf.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.processors import apply_processors, to_list
from protrend.utils.miscellaneous import is_null


class SourceTransformer(CollectfTransformer):
    name = 'collectf'
    type = 'database'
    url = 'http://collectf.org/'
    doi = '10.1093/nar/gkt1123'
    authors = ['Sefa Kiliç', 'Elliot R White', 'Dinara M Sagitova', 'Joseph P Cornish', 'Ivan Erill']
    description = 'CollecTF: a database of experimentally validated transcription factor-binding sites in Bacteria'

    default_node = Source
    default_node_factors = ('name',)
    default_order = 100
    columns = {'protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'}

    def transform(self):
        collectf = dict(name=[self.name],
                        type=[self.type],
                        url=[self.url],
                        doi=[self.doi],
                        authors=[self.authors],
                        description=[self.description])

        df = pd.DataFrame(collectf, index=[0])

        self._stack_transformed_nodes(df)

        return df


class OrganismToSourceConnector(CollectfConnector):
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

        # http://www.collectf.org/browse/view_motif_reports_by_taxonomy/
        url = []
        external_identifier = []
        for tax_id in organism['taxonomy']:
            if not is_null(tax_id):
                url.append(f'http://www.collectf.org/browse/view_motif_reports_by_taxonomy/{tax_id}')
                external_identifier.append(tax_id)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['view_motif_reports_by_taxonomy'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryFamilyToSourceConnector(CollectfConnector):
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

        kwargs = dict(url=['http://collectf.org/browse/browse_by_TF/'] * size,
                      external_identifier=[None] * size,
                      key=['browse_by_TF'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToSourceConnector(CollectfConnector):
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

        # http://www.collectf.org/uniprot/
        url = []
        external_identifier = []
        for uniprot_accession in regulator['uniprot_accession']:
            if not is_null(uniprot_accession):
                url.append(f'http://www.collectf.org/uniprot/{uniprot_accession}')
                external_identifier.append(uniprot_accession)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['uniprot'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class OperonToSourceConnector(CollectfConnector):
    default_from_node = Operon
    default_to_node = Source
    default_connect_stack = {'operon': 'integrated_operon.json', 'source': 'integrated_source.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, regulon=to_list)
        operon = operon.explode('regulon')
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://www.collectf.org/uniprot/
        url = []
        external_identifier = []
        for uniprot_accession in operon['regulon']:
            if not is_null(uniprot_accession):
                url.append(f'http://www.collectf.org/uniprot/{uniprot_accession}')
                external_identifier.append(uniprot_accession)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['uniprot'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class GeneToSourceConnector(CollectfConnector):
    default_from_node = Gene
    default_to_node = Source
    default_connect_stack = {'gene': 'integrated_gene.json', 'source': 'integrated_source.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = apply_processors(gene, regulon=to_list)
        gene = gene.explode('regulon')
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://www.collectf.org/uniprot/
        url = []
        external_identifier = []
        for uniprot_accession in gene['regulon']:
            if not is_null(uniprot_accession):
                url.append(f'http://www.collectf.org/uniprot/{uniprot_accession}')
                external_identifier.append(uniprot_accession)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['uniprot'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class TFBSToSourceConnector(CollectfConnector):
    default_from_node = TFBS
    default_to_node = Source
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'source': 'integrated_source.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = apply_processors(tfbs, regulon=to_list)
        tfbs = tfbs.explode('regulon')
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://www.collectf.org/uniprot/
        url = []
        external_identifier = []
        for uniprot_accession in tfbs['regulon']:
            if not is_null(uniprot_accession):
                url.append(f'http://www.collectf.org/uniprot/{uniprot_accession}')
                external_identifier.append(uniprot_accession)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['uniprot'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToSourceConnector(CollectfConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = apply_processors(rin, regulon=to_list)
        rin = rin.explode('regulon')
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        # http://www.collectf.org/uniprot/
        url = []
        external_identifier = []
        for uniprot_accession in rin['regulon']:
            if not is_null(uniprot_accession):
                url.append(f'http://www.collectf.org/uniprot/{uniprot_accession}')
                external_identifier.append(uniprot_accession)
            else:
                url.append(None)
                external_identifier.append(None)

        kwargs = dict(url=url,
                      external_identifier=external_identifier,
                      key=['uniprot'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
