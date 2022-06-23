from Bio import motifs
from Bio.Seq import Seq

def PWM (aligned_seqs):
    seqs = []
    for seq in aligned_seqs:
        seqs.append(Seq(seq))
    m = motifs.create(seqs)
    descriptor = m.counts.normalize(pseudocounts=0.5)
    return descriptor