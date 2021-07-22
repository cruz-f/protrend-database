from typing import List, Union, Type

import pandas as pd

from protrend.bioapis.gene import NCBIGene
from protrend.bioapis.protrein import UniProtProtein, NCBIProtein
from protrend.bioapis.uniprot import map_uniprot_identifiers
from protrend.transform.annotation.dto import GeneDTO
from protrend.utils.miscellaneous import args_length, scale_arg


def fetch_genes(identifiers: List[str],
                cls: Union[Type[NCBIGene], Type[NCBIProtein], Type[UniProtProtein]],
                taxa: List[str],
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


def _find_in_mapping(acc: str, mapping: pd.DataFrame):
    acc_mask = mapping.loc[:, 'From'] == acc
    acc_series: pd.Series = mapping.loc[acc_mask, 'To']

    return list(acc_series)


def _annotate_uniprot(uniprot_protein: UniProtProtein, gene_dto: GeneDTO):
    if uniprot_protein.identifier:
        gene_dto.uniprot_accession.add(uniprot_protein.identifier)
        gene_dto.locus_tag.add(uniprot_protein.locus_tag)
        gene_dto.name.add(uniprot_protein.name)
        gene_dto.synonyms.update(uniprot_protein.synonyms)
        gene_dto.function.add(uniprot_protein.function)
        gene_dto.description.add(uniprot_protein.description)
        gene_dto.sequence.add(uniprot_protein.sequence)

        return 1, uniprot_protein.identifier

    return 0, ''


def _annotate_ncbi_protein(ncbi_protein: NCBIProtein, gene_dto: GeneDTO):
    if ncbi_protein.identifier:

        gene_dto.ncbi_protein.add(ncbi_protein.identifier)
        gene_dto.locus_tag.add(ncbi_protein.locus_tag)
        gene_dto.name.add(ncbi_protein.name)
        gene_dto.synonyms.update(ncbi_protein.synonyms)
        gene_dto.sequence.add(ncbi_protein.sequence)

        refseq = ''
        genbank = ''

        if ncbi_protein.is_refseq():
            gene_dto.refseq_accession.add(ncbi_protein.refseq_accession)
            refseq = ncbi_protein.refseq_accession
        else:
            gene_dto.genbank_accession.add(ncbi_protein.genbank_accession)
            genbank = ncbi_protein.genbank_accession

        return 1, ncbi_protein.identifier, refseq, genbank

    return 0, '', '', ''


def _annotate_ncbi_gene(ncbi_gene: NCBIGene, gene_dto: GeneDTO):
    if ncbi_gene.identifier:
        gene_dto.ncbi_gene.add(ncbi_gene.identifier)
        gene_dto.locus_tag.add(ncbi_gene.locus_tag)
        gene_dto.name.add(ncbi_gene.name)
        gene_dto.synonyms.update(ncbi_gene.synonyms)
        gene_dto.function.add(ncbi_gene.function)
        gene_dto.position_left.add(ncbi_gene.position_left)
        gene_dto.position_right.add(ncbi_gene.position_right)

        return 1, ncbi_gene.identifier

    return 0, ''


def _annotate_gene(uniprot_protein: UniProtProtein, gene_dto: GeneDTO):
    gene_dto.locus_tag.add(uniprot_protein.locus_tag)
    gene_dto.name.add(uniprot_protein.name)
    gene_dto.synonyms.update(uniprot_protein.synonyms)


def annotate_genes(loci: List[str] = None,
                   names: List[str] = None,
                   taxa: List[str] = None,
                   uniprot_proteins: List[str] = None,
                   ncbi_proteins: List[str] = None,
                   ncbi_genbanks: List[str] = None,
                   ncbi_refseqs: List[str] = None,
                   ncbi_genes: List[str] = None) -> GeneDTO:
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

    :param loci:
    :param names:
    :param taxa:
    :param uniprot_proteins:
    :param ncbi_proteins:
    :param ncbi_genbanks:
    :param ncbi_refseqs:
    :param ncbi_genes:
    :return:
    """

    size = args_length(loci, names, taxa, uniprot_proteins, ncbi_proteins, ncbi_genbanks, ncbi_refseqs, ncbi_genes)

    loci = scale_arg(loci, size)
    names = scale_arg(names, size)
    taxa = scale_arg(taxa, size)
    uniprot_proteins = scale_arg(uniprot_proteins, size)
    ncbi_proteins = scale_arg(ncbi_proteins, size)
    ncbi_genbanks = scale_arg(ncbi_genbanks, size)
    ncbi_refseqs = scale_arg(ncbi_refseqs, size)
    ncbi_genes = scale_arg(ncbi_genes, size)

    uniprot_proteins = fetch_genes(identifiers=uniprot_proteins,
                                   cls=UniProtProtein,
                                   taxa=taxa,
                                   loci=loci,
                                   names=names)

    if ncbi_genbanks[0]:
        ncbi_proteins = fetch_genes(identifiers=ncbi_genbanks,
                                    cls=NCBIProtein,
                                    taxa=taxa,
                                    loci=loci,
                                    names=names,
                                    is_genbank=True)

    elif ncbi_refseqs[0]:
        ncbi_proteins = fetch_genes(identifiers=ncbi_refseqs,
                                    cls=NCBIProtein,
                                    taxa=taxa,
                                    loci=loci,
                                    names=names,
                                    is_refseq=True)

    else:
        ncbi_proteins = fetch_genes(identifiers=ncbi_proteins,
                                    cls=NCBIProtein,
                                    taxa=taxa,
                                    loci=loci,
                                    names=names)

    ncbi_genes = fetch_genes(identifiers=ncbi_genes,
                             cls=NCBIGene,
                             taxa=taxa,
                             loci=loci,
                             names=names)

    # from acc to ncbi protein
    accessions = [protein.identifier for protein in uniprot_proteins if protein.identifier]
    uniprot_ncbi_proteins = map_uniprot_identifiers(accessions, from_='ACC', to='P_GI')
    uniprot_ncbi_refseqs = map_uniprot_identifiers(accessions, from_='ACC', to='P_REFSEQ_AC')
    uniprot_ncbi_genbanks = map_uniprot_identifiers(accessions, from_='ACC', to='EMBL')

    gene_dtos = []

    for uniprot_protein, ncbi_protein, ncbi_gene in zip(uniprot_proteins, ncbi_proteins, ncbi_genes):

        gene_dto = GeneDTO()

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

        gene_dto.ncbi_protein.update(uniprot_ncbi_protein)
        gene_dto.refseq_accession.update(uniprot_ncbi_refseq)
        gene_dto.genbank_accession.update(uniprot_ncbi_genbank)

        if ncbi_protein_id in uniprot_ncbi_protein:
            score += 1

        if ncbi_protein_ref in uniprot_ncbi_refseq:
            score += 1

        if ncbi_protein_gen in uniprot_ncbi_genbank:
            score += 1

        gene_dto.annotation_score = score

        gene_dtos.append(gene_dto)

    return gene_dtos