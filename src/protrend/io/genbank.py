from collections import defaultdict
from typing import Dict, Union, List

import pandas as pd
from Bio import SeqIO

from protrend.utils import SetList


# ----------------------------------------
# Read GenBank Record
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


def read_genbank(file: str) -> pd.DataFrame:
    sequence = SeqIO.read(file, "genbank")

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
