from abc import abstractmethod
from collections import defaultdict
from typing import Union, Callable, Dict, List

import pandas as pd
from Bio.SeqRecord import SeqRecord

from protrend.annotation import (GeneDTO, annotate_genes,
                                 OrganismDTO, annotate_organisms,
                                 PublicationDTO, annotate_publications,
                                 EffectorDTO, annotate_effectors, PathwayDTO, annotate_pathways)
from protrend.log import ProtrendLogger
from .transformations import drop_duplicates, drop_empty_string
from .transformer import Transformer
from protrend.utils import SetList, apply_processors
from protrend.utils.processors import to_list_nan, protrend_hash


def get_values(df, col):
    series = df.get(col)

    if series is None:
        return

    return series.to_list()


class EffectorMixIn:

    @staticmethod
    def annotate_effectors(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        effectors = [EffectorDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(effectors)} effectors')

        names = get_values(df, 'name')

        ProtrendLogger.log.info('Annotating with the following params: name')

        annotate_effectors(dtos=effectors, names=names)

        effectors_dict = [dto.to_dict() for dto in effectors]

        effectors_df = pd.DataFrame(effectors_dict)
        return effectors_df


class GeneMixIn:

    @staticmethod
    def annotate_genes(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        genes = [GeneDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(genes)} genes')

        loci = get_values(df, 'locus_tag')
        names = get_values(df, 'name')
        taxa = get_values(df, 'ncbi_taxonomy')
        uniprot_proteins = get_values(df, 'uniprot_accession')
        ncbi_proteins = get_values(df, 'ncbi_protein')
        ncbi_genbanks = get_values(df, 'genbank_accession')
        ncbi_refseqs = get_values(df, 'refseq_accession')
        ncbi_genes = get_values(df, 'ncbi_gene')

        iterator = zip(
            ('locus_tag', 'name', 'ncbi_taxonomy', 'uniprot_accession', 'ncbi_protein', 'genbank_accession',
             'refseq_accession', 'ncbi_gene'),
            (loci, names, taxa, uniprot_proteins, ncbi_proteins, ncbi_genbanks, ncbi_refseqs, ncbi_genes)
        )

        params = [param for param, value in iterator if value is not None]
        params = ','.join(params)

        ProtrendLogger.log.info(f'Annotating with the following params: {params}')

        annotate_genes(dtos=genes,
                       loci=loci,
                       names=names,
                       taxa=taxa,
                       uniprot_proteins=uniprot_proteins,
                       ncbi_proteins=ncbi_proteins,
                       ncbi_genbanks=ncbi_genbanks,
                       ncbi_refseqs=ncbi_refseqs,
                       ncbi_genes=ncbi_genes)

        genes_dict = [dto.to_dict() for dto in genes]

        genes_df = pd.DataFrame(genes_dict)

        strand_mask = (genes_df['strand'] != 'reverse') & (genes_df['strand'] != 'forward')
        genes_df.loc[strand_mask, 'strand'] = None

        return genes_df


class OrganismMixIn:
    species = ['']
    strain = ['']
    ncbi_taxonomy = [None]
    refseq_accession = ['']
    refseq_ftp = ['']
    genbank_accession = ['']
    genbank_ftp = ['']
    ncbi_assembly = [None]
    assembly_accession = ['']
    name = ['']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    @staticmethod
    def annotate_organisms(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        organisms = [OrganismDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(organisms)} organisms')

        identifiers = get_values(df, 'ncbi_taxonomy')
        names = get_values(df, 'name')

        iterator = zip(
            ('ncbi_taxonomy', 'name',),
            (identifiers, names)
        )

        params = [param for param, value in iterator if value is not None]
        params = ', '.join(params)

        ProtrendLogger.log.info(f'Annotating with the following params: {params}')

        annotate_organisms(dtos=organisms, identifiers=identifiers, names=names)

        organisms_dict = [dto.to_dict() for dto in organisms]

        organisms_df = pd.DataFrame(organisms_dict)
        return organisms_df

    def transform(self: Union[Transformer, 'OrganismMixIn']):
        org = dict(name=self.name,
                   species=self.species,
                   strain=self.strain,
                   ncbi_taxonomy=self.ncbi_taxonomy,
                   refseq_accession=self.refseq_accession,
                   refseq_ftp=self.refseq_ftp,
                   genbank_accession=self.genbank_accession,
                   genbank_ftp=self.genbank_ftp,
                   ncbi_assembly=self.ncbi_assembly,
                   assembly_accession=self.assembly_accession)

        df = pd.DataFrame(org, index=list(range(len(self.name))))

        self.stack_transformed_nodes(df)

        return df


class PathwayMixIn:

    @staticmethod
    def annotate_pathways(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        pathways = [PathwayDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(pathways)} pathways')

        names = get_values(df, 'name')

        ProtrendLogger.log.info('Annotating with the following params: name')

        annotate_pathways(dtos=pathways, names=names)

        pathways_dict = [dto.to_dict() for dto in pathways]

        pathways_df = pd.DataFrame(pathways_dict)
        return pathways_df


class PublicationMixIn:

    @staticmethod
    def annotate_publications(df: pd.DataFrame) -> pd.DataFrame:

        input_values = get_values(df, 'input_value')

        publications = [PublicationDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(publications)} publications')

        identifiers = get_values(df, 'pmid')

        ProtrendLogger.log.info('Annotating with the following params: pmid')

        annotate_publications(dtos=publications, identifiers=identifiers)

        publications_dict = [dto.to_dict() for dto in publications]

        publications_df = pd.DataFrame(publications_dict)
        return publications_df


class RegulatoryInteractionMixIn:

    @staticmethod
    def interaction_hash(df: pd.DataFrame) -> pd.DataFrame:
        # filter by organism + regulator + gene + tfbs + effector + regulatory effect
        df2 = apply_processors(df,
                               organism=to_list_nan,
                               regulator=to_list_nan,
                               gene=to_list_nan,
                               tfbs=to_list_nan,
                               effector=to_list_nan,
                               regulatory_effect=to_list_nan)

        ri_series_hash = df2['organism'].copy()
        ri_series_hash += df2['regulator'].copy()
        ri_series_hash += df2['gene'].copy()
        ri_series_hash += df2['tfbs'].copy()
        ri_series_hash += df2['effector'].copy()
        ri_series_hash += df2['regulatory_effect'].copy()

        df = df.assign(regulatory_interaction_hash=ri_series_hash)
        df = apply_processors(df, regulatory_interaction_hash=protrend_hash)
        df = drop_duplicates(df=df, subset=['regulatory_interaction_hash'], perfect_match=True)
        df = df.dropna(subset=['regulatory_interaction_hash'])

        return df

    @abstractmethod
    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        pass

    def _transform(self: Union[Transformer, 'RegulatoryInteractionMixIn'],
                   network: pd.DataFrame,
                   organism: pd.DataFrame,
                   organism_key: str,
                   regulator: pd.DataFrame,
                   regulator_key: str,
                   gene: pd.DataFrame,
                   gene_key: str,
                   tfbs: pd.DataFrame = None,
                   tfbs_key: str = None,
                   effector: pd.DataFrame = None,
                   effector_key: str = None,
                   regulatory_effect_processor: Callable = None) -> pd.DataFrame:
        network = self.transform_network(network)
        organism = self.transform_organism(organism)
        regulator = self.transform_regulator(regulator)
        gene = self.transform_gene(gene)
        if tfbs is not None:
            tfbs = self.transform_tfbs(tfbs)
        if effector is not None:
            effector = self.transform_effector(effector)

        regulatory_interaction = pd.merge(network, organism, on=organism_key)
        regulatory_interaction = pd.merge(regulatory_interaction, regulator, on=regulator_key)
        regulatory_interaction = pd.merge(regulatory_interaction, gene, on=gene_key)

        if tfbs is not None:
            regulatory_interaction = pd.merge(regulatory_interaction, tfbs, how='left', on=tfbs_key)

        else:
            regulatory_interaction = regulatory_interaction.assign(tfbs=None)

        if effector is not None:
            regulatory_interaction = pd.merge(regulatory_interaction, effector, how='left', on=effector_key)

        else:
            regulatory_interaction = regulatory_interaction.assign(effector=None)

        if regulatory_effect_processor is not None:
            regulatory_interaction = apply_processors(regulatory_interaction,
                                                      regulatory_effect=regulatory_effect_processor)

        regulatory_interaction = self.interaction_hash(regulatory_interaction)
        return regulatory_interaction


class SourceMixIn:
    name = ['']
    type = ['']
    url = ['']
    doi = ['']
    authors = [[]]
    description = ['']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self: Union[Transformer, 'SourceMixIn']):
        db = dict(name=self.name,
                  type=self.type,
                  url=self.url,
                  doi=self.doi,
                  authors=self.authors,
                  description=self.description)

        df = pd.DataFrame(db, index=list(range(len(self.name))))

        self.stack_transformed_nodes(df)

        return df


class TFBSMixIn:

    @staticmethod
    def site_hash(df: pd.DataFrame) -> pd.DataFrame:
        # filter by organism + sequence + strand + start + stop + length
        df2 = apply_processors(df,
                               organism=to_list_nan,
                               sequence=to_list_nan,
                               strand=to_list_nan,
                               start=to_list_nan,
                               stop=to_list_nan,
                               length=to_list_nan)

        ri_series_hash = df2['organism'].copy()
        ri_series_hash += df2['sequence'].copy()
        ri_series_hash += df2['strand'].copy()
        ri_series_hash += df2['start'].copy()
        ri_series_hash += df2['stop'].copy()
        ri_series_hash += df2['length'].copy()

        df = df.assign(site_hash=ri_series_hash)
        df = apply_processors(df, site_hash=protrend_hash)
        df = drop_duplicates(df=df, subset=['site_hash'])
        df = df.dropna(subset=['site_hash'])
        df = drop_empty_string(df, 'site_hash')

        return df


# ----------------------------------------
# Transform GenBank Record
# ----------------------------------------
def get_gene_name(qualifiers: Dict) -> Union[str, None]:
    gene: List[str] = qualifiers.get('gene', [None])

    if gene[0]:
        return gene[0].lower()

    return


def get_locus(qualifiers: Dict) -> Union[str, None]:
    locus: List[str] = qualifiers.get('locus_tag', [None])

    if locus[0]:
        return locus[0].lower()

    return


def get_genbank(qualifiers: Dict) -> Union[str, None]:
    gb: List[str] = qualifiers.get('protein_id', [None])

    if gb[0]:
        return gb[0]

    return


def get_uniprot(qualifiers: Dict) -> Union[str, None]:
    xrefs: List[str] = qualifiers.get('db_xref', [])

    for xref in xrefs:

        if 'UniProtKB/Swiss-Prot:' in xref:
            return xref.replace('UniProtKB/Swiss-Prot:', '').rstrip().lstrip()

    return


class SequenceMixIn:

    @staticmethod
    def transform_sequence(sequence: SeqRecord) -> pd.DataFrame:

        sequence_idx = defaultdict(SetList)

        for feature in sequence.features:

            if feature.type == 'CDS':

                gene_name = get_gene_name(feature.qualifiers)

                if gene_name:
                    sequence_idx[gene_name].append(gene_name)

                    locus = get_locus(feature.qualifiers)
                    sequence_idx[gene_name].append(locus)

                    gb_acc = get_genbank(feature.qualifiers)
                    sequence_idx[gene_name].append(gb_acc)

                    uniprot_acc = get_uniprot(feature.qualifiers)
                    sequence_idx[gene_name].append(uniprot_acc)

        # filter out duplicated gene names
        sequence_idx = {i: val for i, val in enumerate(sequence_idx.values())
                        if len(val) == 4}

        return pd.DataFrame.from_dict(sequence_idx,
                                      orient='index',
                                      columns=['name_lower', 'locus_tag',
                                               'genbank_accession', 'uniprot_accession'])
