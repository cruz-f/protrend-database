from typing import List

from Bio import motifs
from Bio.Seq import Seq


def consensus_sequence(aligned_sequences: List[str]) -> str:
    sequences = [Seq(sequence) for sequence in aligned_sequences]
    motif = motifs.create(sequences, alphabet='ATGC-')
    return str(motif.consensus)
