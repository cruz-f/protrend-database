import re
from typing import List, Union, Type, TYPE_CHECKING

import pandas as pd
from tqdm import tqdm

from protrend.bioapis import NCBIGene, UniProtProtein, NCBIProtein, map_uniprot_identifiers
from protrend.log import ProtrendLogger
from protrend.utils.miscellaneous import args_length, scale_arg, is_null
from ..utils import SetList

if TYPE_CHECKING:
    from .dto import GeneDTO


locus_tag_pattern = re.compile('^[A-Za-z]+[0-9]{3,}$')


def _map_accession(acc: str, mapping: pd.DataFrame) -> List[str]:
    mask = mapping['from'] == acc
    return mapping.loc[mask, 'to'].to_list()


def _fetch_genes(identifiers: List[str],
                 cls: Union[Type[NCBIGene], Type[NCBIProtein], Type[UniProtProtein]],
                 taxa: List[str],
                 loci: List[str],
                 names: List[str],
                 is_refseq: bool = False,
                 is_genbank: bool = False) -> List[Union[NCBIGene, NCBIProtein, UniProtProtein]]:
    if not identifiers:
        return []

    genes = []

    for identifier, taxonomy, locus_tag, name in tqdm(zip(identifiers, taxa, loci, names),
                                                      desc='gene',
                                                      total=len(taxa)):
        if not is_null(identifier):

            if is_refseq:
                gene = cls(refseq_accession=identifier)

            elif is_genbank:
                gene = cls(genbank_accession=identifier)

            else:
                gene = cls(identifier=identifier)

        else:
            gene = cls(taxonomy=taxonomy, locus_tag=locus_tag, name=name)

        gene.fetch()
        genes.append(gene)

    return genes


def _annotate_uniprot(uniprot_protein: UniProtProtein, gene_dto: 'GeneDTO'):
    if uniprot_protein.identifier:
        gene_dto.uniprot_accession.append(uniprot_protein.identifier)
        gene_dto.locus_tag.append(uniprot_protein.locus_tag)
        gene_dto.name.append(uniprot_protein.name)
        gene_dto.synonyms.extend(uniprot_protein.synonyms)
        gene_dto.function.append(uniprot_protein.function)
        gene_dto.description.append(uniprot_protein.description)
        gene_dto.protein_sequence.append(uniprot_protein.sequence)

        return uniprot_protein.identifier

    return ''


def _annotate_ncbi_protein(ncbi_protein: NCBIProtein, gene_dto: 'GeneDTO'):
    if ncbi_protein.identifier:

        gene_dto.ncbi_protein.append(ncbi_protein.identifier)
        gene_dto.locus_tag.append(ncbi_protein.locus_tag)
        gene_dto.synonyms.extend(ncbi_protein.synonyms)
        gene_dto.protein_sequence.append(ncbi_protein.sequence)

        if ncbi_protein.is_refseq():
            gene_dto.refseq_accession.append(ncbi_protein.refseq_accession)
        else:
            gene_dto.genbank_accession.append(ncbi_protein.genbank_accession)


def _annotate_ncbi_gene(ncbi_gene: NCBIGene, gene_dto: 'GeneDTO'):
    if ncbi_gene.identifier:
        gene_dto.ncbi_gene.append(ncbi_gene.identifier)
        gene_dto.locus_tag.append(ncbi_gene.locus_tag)

        # noinspection PyUnresolvedReferences
        if ncbi_gene.name != ncbi_gene.locus_tag and gene_dto.name.take_first() == gene_dto.locus_tag.take_first():
            gene_dto.name.insert(0, ncbi_gene.name)

        else:
            gene_dto.name.append(ncbi_gene.name)

        gene_dto.synonyms.extend(ncbi_gene.synonyms)
        gene_dto.function.append(ncbi_gene.function)
        gene_dto.start.append(ncbi_gene.start)
        gene_dto.stop.append(ncbi_gene.stop)
        gene_dto.strand.append(ncbi_gene.strand)

        return ncbi_gene.identifier

    return ''


def _annotate_gene(uniprot_protein: UniProtProtein, gene_dto: 'GeneDTO'):
    gene_dto.locus_tag.append(uniprot_protein.locus_tag)
    gene_dto.name.append(uniprot_protein.name)
    gene_dto.synonyms.extend(uniprot_protein.synonyms)


def _annotate_locus_tag(gene_dto: 'GeneDTO'):
    locus_tag = gene_dto.locus_tag.take_first()
    loci = gene_dto.locus_tag.take_all()
    synonyms = gene_dto.synonyms.take_all()

    loci = loci + synonyms
    loci = set(loci)

    if len(locus_tag) < 5:

        for locus in loci:
            if locus_tag_pattern.match(locus):
                gene_dto.locus_tag = SetList([locus], output='take_first')

                return

        gene_dto.locus_tag = SetList([], output='take_first')

    elif '(deprecated)' in locus_tag:
        locus_tag = locus_tag.replace('(deprecated)', '').rstrip().lstrip()
        gene_dto.locus_tag = SetList([locus_tag], output='take_first')

    return


def annotate_genes(dtos: List['GeneDTO'],
                   loci: List[str] = None,
                   names: List[str] = None,
                   taxa: List[str] = None,
                   uniprot_proteins: List[str] = None,
                   ncbi_proteins: List[str] = None,
                   ncbi_genbanks: List[str] = None,
                   ncbi_refseqs: List[str] = None,
                   ncbi_genes: List[str] = None) -> List['GeneDTO']:
    """
    A common method to annotate a given gene with relevant information from UniProt or NCBI.

    A given gene is annotated as follows:

        - 1º Step:
            - query genes to UniProt by accession, locus tag + taxonomy, locus tag, gene name + taxonomy
            - retrieve UniProt identifiers
            - map all genes to NCBI protein, GenBank and RefSeq identifiers

        - 2º Step:
            - query genes to NCBI protein by protein id, locus tag + taxonomy, locus tag, gene name + taxonomy
            - retrieve NCBI protein identifiers
            - map all genes to UniProt, GenBank and RefSeq identifiers

        - 3º Step:
            - query genes to NCBI gene by gene id locus tag + taxonomy, locus tag, gene name + taxonomy
            - retrieve NCBI gene identifiers
            - map all genes to NCBI protein, GenBank and RefSeq identifiers

        - 4º Step:
            - merge data into a 'GeneDTO'


    :type dtos: List['GeneDTO']
    :type loci: List[str]
    :type names: List[str]
    :type taxa: List[str]
    :type uniprot_proteins: List[str]
    :type ncbi_proteins: List[str]
    :type ncbi_genbanks: List[str]
    :type ncbi_refseqs: List[str]
    :type ncbi_genes: List[str]

    :rtype: List['GeneDTO']
    
    :param dtos: list of 'GeneDTO'. Each gene DTO will be annotated with information from NCBI and UniProt
    :param loci: list of locus tag identifiers to assist with the annotation
    :param names: list of names to assist with the annotation
    :param taxa: list of taxonomy identifiers to assist with the annotation
    :param uniprot_proteins: list of UniProt accessions to assist with the annotation
    :param ncbi_proteins: list of NCBI Protein identifiers to assist with the annotation
    :param ncbi_genbanks: list of NCBI GenBank accessions to assist with the annotation
    :param ncbi_refseqs: list of NCBI RefSeq accessions to assist with the annotation
    :param ncbi_genes: list of NCBI Gene identifiers to assist with the annotation

    :return: list of annotated 'GeneDTO'. This function returns the same list object for convenience
    """
    if not dtos:
        return dtos

    dtos_size = len(dtos)

    size = args_length(loci, names, taxa, uniprot_proteins, ncbi_proteins, ncbi_genbanks, ncbi_refseqs, ncbi_genes)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    loci = scale_arg(loci, size)
    names = scale_arg(names, size)
    taxa = scale_arg(taxa, size)
    uniprot_proteins = scale_arg(uniprot_proteins, size)
    ncbi_proteins = scale_arg(ncbi_proteins, size)
    ncbi_genbanks = scale_arg(ncbi_genbanks, size)
    ncbi_refseqs = scale_arg(ncbi_refseqs, size)
    ncbi_genes = scale_arg(ncbi_genes, size)

    ProtrendLogger.log.info(f'Starting fetch {len(uniprot_proteins)} genes to {UniProtProtein.__name__}')
    uniprot_proteins = _fetch_genes(identifiers=uniprot_proteins,
                                    cls=UniProtProtein,
                                    taxa=taxa,
                                    loci=loci,
                                    names=names)

    if ncbi_genbanks[0]:
        ProtrendLogger.log.info(f'Starting fetch {len(uniprot_proteins)} genes to {NCBIProtein.__name__} using GenBank')
        ncbi_proteins = _fetch_genes(identifiers=ncbi_genbanks,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names,
                                     is_genbank=True)

    elif ncbi_refseqs[0]:
        ProtrendLogger.log.info(f'Starting fetch {len(uniprot_proteins)} genes to {NCBIProtein.__name__} using RefSeq')
        ncbi_proteins = _fetch_genes(identifiers=ncbi_refseqs,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names,
                                     is_refseq=True)

    else:
        ProtrendLogger.log.info(f'Starting fetch {len(uniprot_proteins)} genes to {NCBIProtein.__name__}')
        ncbi_proteins = _fetch_genes(identifiers=ncbi_proteins,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names)

    ProtrendLogger.log.info(f'Starting fetch {len(uniprot_proteins)} genes to {NCBIGene.__name__}')
    ncbi_genes = _fetch_genes(identifiers=ncbi_genes,
                              cls=NCBIGene,
                              taxa=taxa,
                              loci=loci,
                              names=names)

    ProtrendLogger.log.info(f'Finishing fetch genes')

    # using the uniprot mapping identifier tool
    accessions = [protein.identifier for protein in uniprot_proteins if protein.identifier]

    ProtrendLogger.log.info(f'Starting map genes with '
                            f'{len(accessions)} uniprot accessions to ncbi proteins (GI_number)')
    uniprot_ncbi_proteins = map_uniprot_identifiers(accessions, from_db='UniProtKB_AC-ID', to_db='GI_number')

    ProtrendLogger.log.info(f'Starting map genes with '
                            f'{len(accessions)} uniprot accessions to ncbi refseqs (RefSeq_Protein)')
    uniprot_ncbi_refseqs = map_uniprot_identifiers(accessions, from_db='UniProtKB_AC-ID', to_db='RefSeq_Protein')

    ProtrendLogger.log.info(f'Starting map genes with '
                            f'{len(accessions)} uniprot accessions to ncbi genbanks (RefSeq_Nucleotide)')
    uniprot_ncbi_genbanks = map_uniprot_identifiers(accessions, from_db='UniProtKB_AC-ID', to_db='RefSeq_Nucleotide')

    ProtrendLogger.log.info(f'Starting map genes with '
                            f'{len(accessions)} uniprot accessions to ncbi gene (GeneID)')
    uniprot_ncbi_genes = map_uniprot_identifiers(accessions, from_db='UniProtKB_AC-ID', to_db='GeneID')

    for gene_dto, uniprot_protein, ncbi_protein, ncbi_gene in zip(dtos, uniprot_proteins, ncbi_proteins, ncbi_genes):
        uniprot_id = _annotate_uniprot(uniprot_protein, gene_dto)
        _annotate_ncbi_protein(ncbi_protein, gene_dto)
        _annotate_ncbi_gene(ncbi_gene, gene_dto)

        # base annotation, in case none of the annotations work
        _annotate_gene(uniprot_protein, gene_dto)

        # curating the locus tag annotation due to mis annotations of some locus tag
        _annotate_locus_tag(gene_dto)

        uniprot_ncbi_protein = _map_accession(uniprot_id, uniprot_ncbi_proteins)
        uniprot_ncbi_refseq = _map_accession(uniprot_id, uniprot_ncbi_refseqs)
        uniprot_ncbi_genbank = _map_accession(uniprot_id, uniprot_ncbi_genbanks)
        uniprot_ncbi_gene = _map_accession(uniprot_id, uniprot_ncbi_genes)

        gene_dto.ncbi_protein.extend(uniprot_ncbi_protein)
        gene_dto.refseq_accession.extend(uniprot_ncbi_refseq)
        gene_dto.genbank_accession.extend(uniprot_ncbi_genbank)
        gene_dto.ncbi_gene.extend(uniprot_ncbi_gene)

    return dtos
