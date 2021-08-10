from typing import List, Union, Tuple

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source, Organism, Effector, Pathway, RegulatoryFamily, Publication
from protrend.transform.annotation.gene import annotate_genes
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise.settings import RegulatorSettings
from protrend.transform.transformer import Transformer


class RegulatorTransformer(Transformer):

    def __init__(self, settings: RegulatorSettings = None):

        if not settings:
            settings = RegulatorSettings()

        super().__init__(settings)

    def _read_regulon(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('regulon')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam',
                                       'biological_process', 'regulation_effector', 'regulation_regulog',
                                       'regulog', 'taxonomy', 'rna_family', 'effector', 'pathway', 'operon',
                                       'tfbs', 'gene', 'regulator_locus_tag', 'regulator_family',
                                       'regulation_mode', 'transcription_factor', 'tf_family'])

        return df

    def _read_organism(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('organism')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=['protrend_id',
                                       'genome_id', 'name', 'taxonomy', 'url', 'regulon',
                                       'species', 'strain', 'family', 'phylum',
                                       'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                                       'genbank_accession', 'genbank_ftp',
                                       'ncbi_assembly', 'assembly_accession'])

        df = df.rename(columns={'protrend_id': 'organism_protrend_id'})

        df = df[['organism_protrend_id', 'genome_id', 'ncbi_taxonomy']]

        return df

    def _transform_tf(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:

        # filter tfs only
        mask = regulon['regulator_locus_tag'].notnull()
        regulon = regulon[mask]

        regulon = self.drop_duplicates(df=regulon,
                                       subset=['regulator_locus_tag'],
                                       perfect_match=True,
                                       preserve_nan=False)

        apply_processors(rstrip,
                         lstrip,
                         df=regulon,
                         col='regulator_locus_tag')

        apply_processors(rstrip,
                         lstrip,
                         df=regulon,
                         col='name')

        apply_processors(rstrip,
                         lstrip,
                         df=organism,
                         col='genome_id')

        regulon['mechanism'] = ['transcription factor'] * regulon.shape[0]

        df = pd.merge(regulon, organism, left_on='genome', right_on='genome_id')

        df['input_value'] = df['regulator_locus_tag']

        return df

    def _transform_rna(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:

        # filter rna only
        mask = regulon['rfam'].notnull()
        regulon = regulon[mask]

        regulon = self.drop_duplicates(df=regulon,
                                       subset=['rfam', 'genome'],
                                       perfect_match=True,
                                       preserve_nan=False)

        apply_processors(rstrip,
                         lstrip,
                         df=regulon,
                         col='rfam')

        apply_processors(rstrip,
                         lstrip,
                         df=regulon,
                         col='name')

        apply_processors(rstrip,
                         lstrip,
                         df=organism,
                         col='genome_id')

        regulon['mechanism'] = ['small RNA (sRNA)'] * regulon.shape[0]

        df = pd.merge(regulon, organism, left_on='genome', right_on='genome_id')

        df['input_value'] = df['name']

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
        regulon = self._read_regulon()
        organism = self._read_organism()

        # ------------------ regulon of type TF --------------------------------
        tf = self._transform_tf(regulon=regulon, organism=organism)

        loci = tf['regulator_locus_tag'].tolist()
        names = tf['name'].tolist()
        taxa = tf['ncbi_taxonomy'].tolist()

        tf_genes = self._annotate_genes(loci, names, taxa)

        tf_df = pd.merge(tf_genes, tf, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        tf_df['name'] = tf_df['name_annotation']

        tf_df = tf_df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        # ------------------ regulon of type RNA --------------------------------
        rna = self._transform_rna(regulon=regulon, organism=organism)

        loci = [None] * rna.shape[0]
        names = rna['name'].tolist()
        taxa = rna['ncbi_taxonomy'].tolist()

        rna_genes = self._annotate_genes(loci, names, taxa)

        rna_df = pd.merge(rna_genes, rna, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        tf_df['name'] = tf_df['name_annotation']

        tf_df = tf_df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        # --------------------- concat DFs --------------------------------------

        df = pd.concat([tf_df, rna_df], axis=0)

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
        from_identifiers = from_df['protrend_id'].tolist()

        size = len(from_identifiers)

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]
        to_identifiers = [regprecise_id] * size

        kwargs = dict(url=from_df['url'].tolist(),
                      external_identifier=from_df['regulon_id'].tolist(),
                      key=['regulon_id'] * size)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    kwargs=kwargs)

    def _connect_to_organism(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')

        if not from_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        from_identifiers = from_df['protrend_id'].tolist()
        to_identifiers = from_df['organism_protrend_id'].tolist()

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Organism,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def _connect_to_effector(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_effector')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        to_df = read_csv(to_path)

        from_identifiers = []
        to_identifiers = []

        for i, effector_row in enumerate(from_df['effector']):

            if not effector_row:
                continue

            for effector_id in effector_row:
                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['effector_id'].values == effector_id.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Effector,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def _connect_to_pathway(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_pathway')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        to_df = read_csv(to_path)

        from_identifiers = []
        to_identifiers = []

        for i, pathway_row in enumerate(from_df['pathway']):

            if not pathway_row:
                continue

            for effector_id in pathway_row:
                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['pathway_id'].values == effector_id.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Pathway,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def _connect_to_regulatory_family_publication(self) -> Tuple[pd.DataFrame, pd.DataFrame]:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_regulatory_family')

        if not from_path:
            return pd.DataFrame(), pd.DataFrame()

        if not to_path:
            return pd.DataFrame(), pd.DataFrame()

        from_df = read_csv(from_path)
        to_df = read_csv(to_path)

        from_identifiers = []
        to_identifiers = []

        for i, (tf_family, tf, rna_family) in enumerate(zip(from_df['rna_family'],
                                                            from_df['transcription_factor'],
                                                            from_df['tf_family'])):

            if tf_family:
                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['tffamily_id'].values == tf_family.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

            elif tf:
                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['collection_id'].values == tf.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

            elif rna_family:
                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['riboswitch_id'].values == rna_family.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

        size = len(from_identifiers)

        from_reg_to_fam = self.make_connection(size=size,
                                               from_node=self.node,
                                               to_node=RegulatoryFamily,
                                               from_identifiers=from_identifiers,
                                               to_identifiers=to_identifiers)

        to_path_pub = self._connect_stack.get('to_publication')

        if not to_path_pub:
            return from_reg_to_fam, pd.DataFrame()

        from_pub_ids = []
        to_pub_ids = []

        for from_id, to_id in zip(from_identifiers, to_identifiers):

            mask = to_df['protrend_id'] == to_id
            pubmed_row = to_df.loc[mask, 'pubmed']

            if not pubmed_row:
                continue

            for pmid in pubmed_row:
                from_pub_ids.append(from_id)

                pmid_mask = to_path_pub['pmid'].values == pmid.replace(' ', '')
                pub_id = to_path_pub.loc[pmid_mask, 'protrend_id'].iloc[0]

                to_pub_ids.append(pub_id)

        size = len(from_pub_ids)

        from_reg_to_pub = self.make_connection(size=size,
                                               from_node=self.node,
                                               to_node=Publication,
                                               from_identifiers=from_pub_ids,
                                               to_identifiers=to_pub_ids)

        return from_reg_to_fam, from_reg_to_pub

    def connect(self):

        connection = self._connect_to_source()
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, connection)

        connection = self._connect_to_organism()
        df_name = f'connected_{self.node.node_name()}_{Organism.node_name()}'
        self.stack_csv(df_name, connection)

        connection = self._connect_to_effector()
        df_name = f'connected_{self.node.node_name()}_{Effector.node_name()}'
        self.stack_csv(df_name, connection)

        connection = self._connect_to_pathway()
        df_name = f'connected_{self.node.node_name()}_{Pathway.node_name()}'
        self.stack_csv(df_name, connection)

        fam_connection, pub_connection = self._connect_to_regulatory_family_publication()

        df_name = f'connected_{self.node.node_name()}_{RegulatoryFamily.node_name()}'
        self.stack_csv(df_name, fam_connection)

        df_name = f'connected_{self.node.node_name()}_{Publication.node_name()}'
        self.stack_csv(df_name, pub_connection)
