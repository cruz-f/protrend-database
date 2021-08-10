from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source, Organism
from protrend.transform.annotation.gene import annotate_genes
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, take_first
from protrend.transform.regprecise.settings import GeneSettings
from protrend.transform.transformer import Transformer
from protrend.utils.miscellaneous import take_last, flatten_list


class GeneTransformer(Transformer):

    def __init__(self, settings: GeneSettings = None):

        if not settings:
            settings = GeneSettings()

        super().__init__(settings)

    def _read_gene(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('gene')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'])

        return df

    def _read_regulator(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('regulator')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=['protrend_id',
                                       'organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
                                       'mechanism',
                                       'regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam',
                                       'biological_process', 'regulation_effector', 'regulation_regulog',
                                       'regulog', 'taxonomy', 'rna_family', 'effector', 'pathway', 'operon',
                                       'tfbs', 'gene', 'regulator_locus_tag', 'regulator_family',
                                       'regulation_mode', 'transcription_factor', 'tf_family',
                                       'locus_tag', 'synonyms', 'function', 'description',
                                       'ncbi_gene', 'ncbi_protein',
                                       'genbank_accession', 'refseq_accession', 'uniprot_accession',
                                       'sequence', 'strand', 'position_left', 'position_right', 'annotation_score'])

        df = df.rename(columns={'protrend_id': 'regulator_protrend_id'})

        df = df[['organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
                 'regulator_protrend_id', 'regulon_id']]

        return df

    @staticmethod
    def _transform_gene(gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:

        apply_processors(rstrip,
                         lstrip,
                         df=gene,
                         col='locus_tag')

        apply_processors(rstrip,
                         lstrip,
                         df=gene,
                         col='name')

        aggregation_functions = {'name': take_last,
                                 'function': take_last,
                                 'url': set,
                                 'regulon': flatten_list,
                                 'operon': flatten_list,
                                 'tfbs': flatten_list}

        gene = gene.groupby(gene['locus_tag']).aggregate(aggregation_functions)
        gene = gene.reset_index()

        gene['regulon_id'] = gene['regulon']

        # keeping only one, since we only want to get the ncbi taxonomy of each gene.
        apply_processors(take_first,
                         df=gene,
                         col='regulon_id')

        df = pd.merge(gene, regulator, on='regulon_id')

        df['input_value'] = df['locus_tag']

        return df

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[int]):

        dtos = [GeneDTO(input_value=locus) for locus in loci]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa)

        for dto, name in zip(dtos, names):
            dto.synonyms.append(name)

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])

        # locus_tag: List[str]
        # name: List[str]
        # synonyms: List[str]
        # function: List[str]
        # description: List[str]
        # ncbi_gene: List[str]
        # ncbi_protein: List[str]
        # genbank_accession: List[str]
        # refseq_accession: List[str]
        # uniprot_accession: List[str]
        # sequence: List[str]
        # strand: List[str]
        # position_left: List[int]
        # position_right: List[int]
        # annotation_score: int

        genes_cols = ['input_value',
                      'locus_tag', 'name',
                      'synonyms', 'function',
                      'description'
                      'ncbi_gene', 'ncbi_protein',
                      'genbank_accession', 'refseq_accession', 'uniprot_accession',
                      'sequence',
                      'strand',
                      'position_left',
                      'position_right',
                      'annotation_score']

        if genes.empty:
            genes = pd.DataFrame(columns=genes_cols)

        return genes

    def transform(self):
        gene = self._read_gene()
        regulator = self._read_regulator()

        gene = self._transform_gene(gene=gene, regulator=regulator)

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()

        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        df['locus_tag'] = df['locus_tag_annotation']
        df['name'] = df['name_annotation']
        df['function'] = df['function_annotation']

        df = df.drop(['input_value',
                      'name_annotation', 'name_regprecise',
                      'locus_tag_annotation',
                      'function_annotation', 'function_regprecise'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_source')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)

        from_identifiers = []
        kwargs = defaultdict(list)

        for _, row in from_df.iterrows():

            from_id = row['protrend_id']
            urls = row['url']
            regulons = row['regulon']

            for url, regulon in zip(urls, regulons):
                from_identifiers.append(from_id)
                kwargs['url'].append(url)
                kwargs['external_identifier'].append(regulon)
                kwargs['key'].append('regulon_id')

        size = len(from_identifiers)

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]
        to_identifiers = [regprecise_id] * size

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    kwargs=kwargs)

    def _connect_to_organism(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_organism')

        if not from_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        to_df = read_csv(to_path)

        from_identifiers = []
        to_identifiers = []

        for i, ncbi in enumerate(from_df['ncbi_taxonomy']):

            if not ncbi:
                continue

            from_id = from_df['protrend_id'].iloc[i]
            from_identifiers.append(from_id)

            mask = to_df['ncbi_taxonomy'].values == ncbi.replace(' ', '')
            to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
            to_identifiers.append(to_id)

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Organism,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def connect(self):

        connection = self._connect_to_source()
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, connection)

        connection = self._connect_to_organism()
        df_name = f'connected_{self.node.node_name()}_{Organism.node_name()}'
        self.stack_csv(df_name, connection)
