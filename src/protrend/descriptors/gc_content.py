from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqUtils import GC

def GCContent (aligned_seqs):
    seqs = []
    for seq in aligned_seqs:
        seqs.append(Seq(seq))
    m = motifs.create(seqs, alphabet='ATGC-')
    descriptor = GC(m.consensus)
    return descriptor