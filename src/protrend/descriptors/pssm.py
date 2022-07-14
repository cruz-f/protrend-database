from Bio import motifs
from Bio.Seq import Seq

def PSSM (aligned_seqs):
    seqs = []
    for seq in aligned_seqs:
        seqs.append(Seq(seq))
    m = motifs.create(seqs, alphabet='ATGC-')
    pwm = m.counts.normalize(pseudocounts=0.5)
    descriptor = pwm.log_odds()
    return descriptor