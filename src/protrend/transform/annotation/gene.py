from typing import List, Union, Type

import pandas as pd

from protrend.bioapis.gene import NCBIGene
from protrend.bioapis.protrein import UniProtProtein, NCBIProtein
from protrend.bioapis.uniprot import map_uniprot_identifiers
from protrend.transform.dto import GeneDTO
from protrend.utils.miscellaneous import args_length, scale_arg


def _find_in_mapping(acc: str, mapping: pd.DataFrame):
    acc_mask = mapping.loc[:, 'From'] == acc
    acc_series: pd.Series = mapping.loc[acc_mask, 'To']

    return list(acc_series)


def _fetch_genes(identifiers: List[str],
                 cls: Union[Type[NCBIGene], Type[NCBIProtein], Type[UniProtProtein]],
                 taxa: List[int],
                 loci: List[str],
                 names: List[str],
                 is_refseq: bool = False,
                 is_genbank: bool = False) -> List[Union[NCBIGene, NCBIProtein, UniProtProtein]]:
    genes = []

    if identifiers[0] is None:

        for taxonomy, locus_tag, name in zip(taxa, loci, names):
            gene = cls(taxonomy=taxonomy, locus_tag=locus_tag, name=name)
            gene.fetch()
            genes.append(gene)

    else:

        for identifier in identifiers:

            if is_refseq:
                gene = cls(refseq_accession=identifier)

            elif is_genbank:
                gene = cls(genbank_accession=identifier)

            else:
                gene = cls(identifier=identifier)

            gene.fetch()
            genes.append(gene)

    return genes


def _annotate_uniprot(uniprot_protein: UniProtProtein, gene_dto: GeneDTO):
    if uniprot_protein.identifier:
        gene_dto.uniprot_accession.append(uniprot_protein.identifier)
        gene_dto.locus_tag.append(uniprot_protein.locus_tag)
        gene_dto.name.append(uniprot_protein.name)
        gene_dto.synonyms.extend(uniprot_protein.synonyms)
        gene_dto.function.append(uniprot_protein.function)
        gene_dto.description.append(uniprot_protein.description)
        gene_dto.sequence.append(uniprot_protein.sequence)

        return 1, uniprot_protein.identifier

    return 0, ''


def _annotate_ncbi_protein(ncbi_protein: NCBIProtein, gene_dto: GeneDTO):
    if ncbi_protein.identifier:

        gene_dto.ncbi_protein.append(ncbi_protein.identifier)
        gene_dto.locus_tag.append(ncbi_protein.locus_tag)
        gene_dto.name.append(ncbi_protein.name)
        gene_dto.synonyms.extend(ncbi_protein.synonyms)
        gene_dto.sequence.append(ncbi_protein.sequence)

        refseq = ''
        genbank = ''

        if ncbi_protein.is_refseq():
            gene_dto.refseq_accession.append(ncbi_protein.refseq_accession)
            refseq = ncbi_protein.refseq_accession
        else:
            gene_dto.genbank_accession.append(ncbi_protein.genbank_accession)
            genbank = ncbi_protein.genbank_accession

        return 1, ncbi_protein.identifier, refseq, genbank

    return 0, '', '', ''


def _annotate_ncbi_gene(ncbi_gene: NCBIGene, gene_dto: GeneDTO):
    if ncbi_gene.identifier:
        gene_dto.ncbi_gene.append(ncbi_gene.identifier)
        gene_dto.locus_tag.append(ncbi_gene.locus_tag)
        gene_dto.name.append(ncbi_gene.name)
        gene_dto.synonyms.extend(ncbi_gene.synonyms)
        gene_dto.function.append(ncbi_gene.function)
        gene_dto.position_left.append(ncbi_gene.position_left)
        gene_dto.position_right.append(ncbi_gene.position_right)

        return 1, ncbi_gene.identifier

    return 0, ''


def _annotate_gene(uniprot_protein: UniProtProtein, gene_dto: GeneDTO):
    gene_dto.locus_tag.append(uniprot_protein.locus_tag)
    gene_dto.name.append(uniprot_protein.name)
    gene_dto.synonyms.extend(uniprot_protein.synonyms)


def annotate_genes(dtos: List[GeneDTO],
                   loci: List[str] = None,
                   names: List[str] = None,
                   taxa: List[int] = None,
                   uniprot_proteins: List[str] = None,
                   ncbi_proteins: List[str] = None,
                   ncbi_genbanks: List[str] = None,
                   ncbi_refseqs: List[str] = None,
                   ncbi_genes: List[str] = None) -> List[GeneDTO]:
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
            - merge data into a GeneDTO


    :type dtos: List[GeneDTO]
    :type loci: List[str]
    :type names: List[str]
    :type taxa: List[int]
    :type uniprot_proteins: List[str]
    :type ncbi_proteins: List[str]
    :type ncbi_genbanks: List[str]
    :type ncbi_refseqs: List[str]
    :type ncbi_genes: List[str]

    :rtype: List[GeneDTO]
    
    :param dtos: list of GeneDTO. Each gene DTO will be annotated with information from NCBI and UniProt
    :param loci: list of locus tag identifiers to assist with the annotation
    :param names: list of names to assist with the annotation
    :param taxa: list of taxonomy identifiers to assist with the annotation
    :param uniprot_proteins: list of UniProt accessions to assist with the annotation
    :param ncbi_proteins: list of NCBI Protein identifiers to assist with the annotation
    :param ncbi_genbanks: list of NCBI GenBank accessions to assist with the annotation
    :param ncbi_refseqs: list of NCBI RefSeq accessions to assist with the annotation
    :param ncbi_genes: list of NCBI Gene identifiers to assist with the annotation

    :return: list of annotated GeneDTO. This function returns the same list object for convenience
    """

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

    uniprot_proteins = _fetch_genes(identifiers=uniprot_proteins,
                                    cls=UniProtProtein,
                                    taxa=taxa,
                                    loci=loci,
                                    names=names)

    if ncbi_genbanks[0]:
        ncbi_proteins = _fetch_genes(identifiers=ncbi_genbanks,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names,
                                     is_genbank=True)

    elif ncbi_refseqs[0]:
        ncbi_proteins = _fetch_genes(identifiers=ncbi_refseqs,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names,
                                     is_refseq=True)

    else:
        ncbi_proteins = _fetch_genes(identifiers=ncbi_proteins,
                                     cls=NCBIProtein,
                                     taxa=taxa,
                                     loci=loci,
                                     names=names)

    ncbi_genes = _fetch_genes(identifiers=ncbi_genes,
                              cls=NCBIGene,
                              taxa=taxa,
                              loci=loci,
                              names=names)

    # from acc to ncbi protein
    accessions = [protein.identifier for protein in uniprot_proteins if protein.identifier]
    uniprot_ncbi_proteins = map_uniprot_identifiers(accessions, from_='ACC', to='P_GI')
    uniprot_ncbi_refseqs = map_uniprot_identifiers(accessions, from_='ACC', to='P_REFSEQ_AC')
    uniprot_ncbi_genbanks = map_uniprot_identifiers(accessions, from_='ACC', to='EMBL')

    for gene_dto, uniprot_protein, ncbi_protein, ncbi_gene in zip(dtos, uniprot_proteins, ncbi_proteins, ncbi_genes):

        uniprot_protein_score, uniprot_id = _annotate_uniprot(uniprot_protein, gene_dto)
        ncbi_protein_score, ncbi_protein_id, ncbi_protein_ref, ncbi_protein_gen = _annotate_ncbi_protein(ncbi_protein,
                                                                                                         gene_dto)
        ncbi_gene_score, ncbi_gene_id = _annotate_ncbi_gene(ncbi_gene, gene_dto)

        score = uniprot_protein_score + ncbi_protein_score + ncbi_gene_score

        # base annotation, in case none of the annotations work
        _annotate_gene(uniprot_protein, gene_dto)

        uniprot_ncbi_protein = _find_in_mapping(uniprot_id, uniprot_ncbi_proteins)
        uniprot_ncbi_refseq = _find_in_mapping(uniprot_id, uniprot_ncbi_refseqs)
        uniprot_ncbi_genbank = _find_in_mapping(uniprot_id, uniprot_ncbi_genbanks)

        gene_dto.ncbi_protein.extend(uniprot_ncbi_protein)
        gene_dto.refseq_accession.extend(uniprot_ncbi_refseq)
        gene_dto.genbank_accession.extend(uniprot_ncbi_genbank)

        if ncbi_protein_id in uniprot_ncbi_protein:
            score += 1

        if ncbi_protein_ref in uniprot_ncbi_refseq:
            score += 1

        if ncbi_protein_gen in uniprot_ncbi_genbank:
            score += 1

        gene_dto.annotation_score = score

    return dtos
